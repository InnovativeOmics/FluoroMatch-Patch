#!/usr/bin/env Rscript
library(data.table)
library(SearchTrees)
library(comprehenr)
library(mzR)
library(msdata)
library(txtplot)

arguments = list(
    #  fn_NegID = "NegIDed_FIN_B_small.csv"
    # fn_NegID = "AIF/3MLW_TargetList.csv"
    # fn_NegID = "NegIDed_FIN_small.csv"
    fn_NegID = "NegIDed_FIN.csv"
    ,fn_mzxml = "FoamPool_Targeted_Neg.mzXML"
    # ,fn_mzxml = "AIF/220506_3MLW_NEW5000X_1_ALL_IONS.mzXML"
    ,path_to_output_folder = "."
    ,fn_csv_output = "EXAMPLE_CSV_OUTPUT.csv" # this is converted from mzxml data
    ,fn_eic_output = "EXAMPLE_EIC_OUTPUT.csv" # this is where the EIC is stored using feature table
    ,mztol = 0.005
    ,rttol = 0.5
    # ,NegID_Cols = c(7,8,13) #mz, rt, rowID
    ,NegID_Cols = c(6,7,12) #mz, rt, rowID
    # ,NegID_Cols = c(1,2,3) #mz, rt, rowID
    # ,NegID_SkipRows = c(1) #c() means don't skip, c(1) means skip the first row
    ,NegID_SkipRows = c() #c() means don't skip, c(1) means skip the first row
    ,Min_EIC_Len = 10
    ,precision_mz = 5
    ,precision_rt = 3
    ,min_intensity = 1000
    ,use_min_i_filter = TRUE
    ,cols = c("rt", "mz", "intensity")
    ,AggFUN = max #min, mean, or max: sum will sum the mz and rt would need more code
    ,UseAgg = TRUE 
)

NegID <- function(fn, columns, skip){
    # Read NegID file with columns
    # Allow for skipping extra info column headers
    data <- fread(fn)
    cols <- data[, ..columns]
    if(length(skip) > 0){ 
        cols = cols[-skip,]
    }
    cols <- as.matrix(cols)
    return(cols)
}

EIC_Tree <- function(tree, trt, rttol, tmz, mztol){
    # Return the indicies within a bounding box (mz, rt)
    v = rectLookup(tree, xlim = c(trt-rttol, trt+rttol)
                       , ylim = c(tmz-mztol, tmz+mztol))
    return(v)
}

getCEsize <- function(AllIons, n){
    # Recursive function
    # Get the number of unique collision energies (including NA)
    # Checks 1, 2, 4, 8, until # of ce is less than n
    ce = unique(header(AllIons,1:n)$collisionEnergy)
    if (length(ce) == n) {
        return(getCEsize(AllIons, n*2))
    } else {
        return(length(ce))
    }

}

getSpectra <- function(arguments, AllIons, PEAKS, data, CEs, rt, i){
    # makes a mini matrix of the spectra
    CEi = CEs[i]
    p = PEAKS[[i]]
    if( arguments$use_min_i_filter ){
        p = p[ arguments$min_intensity < p[,2], ]
    }

    m = matrix(nrow = nrow(p), ncol = length(arguments$cols))
    m[,1] = rt[i]
    m[,2:3] = p
    # If needed, set extra column data
    colnames(m) = arguments$cols
    return(m)
}

getAllSpectras <- function(arguments, AllIons, ce_i, CEsize){
    # creates matrix for all scans for a particular CE
    end = length(AllIons)
    # assumes alternating collision energies: 0,10,35,0,10,35
    # potential source of error
    CEs = seq(ce_i, end, CEsize)
    PEAKS = lapply(seq_along(CEs), function(i) peaks(AllIons, CEs[i]))
    data = header(AllIons,CEs)
    rt = data$retentionTime/60
    # If needed, get extra column data
    
    MS = lapply(seq_along(CEs), function(i) getSpectra(arguments, AllIons, PEAKS, data, CEs, rt, i))
    M = do.call(rbind, MS)
    return(M)
}

getNearestRt <- function(rts, rt){
    # return the nearest retention time
    closest = which.min(abs(rts - rt))
    return(rts[closest])
}

generateEICforNegID <- function(CEs, RTs, trees, df_NegID, NegID_row, mztol, rttol){
    # Generate the EIC for a row of NegID
    tmz = df_NegID[NegID_row,1]
    trt = df_NegID[NegID_row,2]
    trt_CE0 = getNearestRt(RTs[[1]], trt)
    EIC_Indices = EIC_Tree(trees[[1]], trt, rttol, tmz, mztol)
    EIC_negID = CEs[[1]][EIC_Indices,]

    return( EIC_negID )
}

printEICs <- function(arguments, EIC_negID, df_NegID, NegID_row){
    # when Min_EIC_Len is 1, skips printing when there is only one row
    if (length(EIC_negID) > length(arguments$cols) * arguments$Min_EIC_Len){
        # rowsum(EIC_negID, EIC_negID[,1])
        if(arguments$UseAgg){
            EIC_negID = aggregate(EIC_negID, by = list(EIC_negID[,1])
                            , FUN = arguments$AggFUN)[,arguments$cols]
        }
        EIC_negID = cbind(EIC_negID, df_NegID[NegID_row,3])
        EIC_negID = cbind(EIC_negID, arguments$fn_mzxml)
        EIC_negID = cbind(EIC_negID, 1)

        colnames(EIC_negID) = c("RT", "mz", "Intensity", "Feature", "File", "Zoom")
        write.table(EIC_negID[,c("Feature", "RT", "Intensity", "mz", "File", "Zoom")]
            , file = arguments$fn_eic_output_wpath, row.names = FALSE, quote = FALSE
            , col.names = FALSE, append = TRUE, sep = ",")
    }

}

calculateEICs <- function(arguments, CEs, RTs, trees, df_NegID){
    # generate the negID EICs
    # write EIC data to file
    for(NegID_row in 1:nrow(df_NegID)){
        print(NegID_row)
        EIC_negID = generateEICforNegID(CEs, RTs, trees, df_NegID, NegID_row, arguments$mztol, arguments$rttol)
        printEICs(arguments, EIC_negID, df_NegID, NegID_row)
    }
}

printCSVheader <- function(fn_csv_output){
    # Write the header for the csv file
    cat(paste("rt", "mz", "intensity", sep=","), file = fn_csv_output, append = FALSE, sep = "\n")
}

printEICheader <- function(fn_eic_output){
    # Write the header for the eic file
    cat(paste('Feature', 'RT', 'Intensity', 'mz', 'File', 'Zoom', sep=","), file = fn_eic_output, append = FALSE, sep = "\n")
}

extractCE2 <- function( arguments ){
    # Convert an mzxml file with multiple collision energies to ms2 with the fragments
    arguments$fn_mzxml_wpath = file.path(arguments$path_to_output_folder, arguments$fn_mzxml)
    arguments$fn_csv_output_wpath = file.path(arguments$path_to_output_folder, arguments$fn_csv_output)
    arguments$fn_eic_output_wpath = file.path(arguments$path_to_output_folder, arguments$fn_eic_output)
    AllIons <- openMSfile(arguments$fn_mzxml_wpath)  
    CEsize = getCEsize(AllIons,1)
    df_NegID <- NegID(arguments$fn_NegID, arguments$NegID_Cols, arguments$NegID_SkipRows)
    class(df_NegID) <- "numeric"

    # printCSVheader(arguments$fn_csv_output_wpath)
    printEICheader(arguments$fn_eic_output_wpath)

    CEs = to_list(for(i in 1:CEsize) getAllSpectras(arguments, AllIons, i, CEsize))
    RTs = to_list(for(i in 1:CEsize) unique(CEs[[i]][,1]))
    trees = to_list(for(i in 1:CEsize) createTree(CEs[[i]], treeType = "quad", dataType = "point", columns = 1:2))

    CEs[[1]][,1] = round(CEs[[1]][,1],arguments$precision_rt)
    CEs[[1]][,2] = round(CEs[[1]][,2],arguments$precision_mz)
    CEs[[1]][,3] = round(CEs[[1]][,3],0)

    # write.table(CEs[[1]][,c("rt", "mz", "intensity")], file = arguments$fn_csv_output_wpath
    #     , row.names = FALSE, quote = FALSE, col.names = FALSE, append = TRUE, sep = ",")
    
    calculateEICs(arguments, CEs, RTs, trees, df_NegID)
}

extractCE2( arguments )