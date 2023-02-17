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
    ,path_to_output_folder = "eic_testing"
    ,fn_mzxml = "FoamPool_Targeted_Neg.mzXML"
    # ,fn_mzxml = "AIF/220506_3MLW_NEW5000X_1_ALL_IONS.mzXML"
    ,fn_MS1_output = "EXAMPLE_MS1_OUTPUT.csv" # this is where the MS1 is stored using feature table
    ,mztol = 0.005
    ,mzZoomWindow = 5
    # ,NegID_Cols = c(7,8,13) #mz, rt, rowID
    ,NegID_Cols = c(6,7,12,4) #mz, rt, rowID, formula
    # ,NegID_Cols = c(1,2,3) #mz, rt, rowID
    # ,NegID_SkipRows = c(1) #c() means don't skip, c(1) means skip the first row
    ,NegID_SkipRows = c() #c() means don't skip, c(1) means skip the first row
    ,precision_mz = 5
    ,precision_rt = 3
    ,min_intensity = 1000
    ,use_min_i_filter = TRUE
    ,cols = c("rt", "mz", "intensity")
    ,IsoThreshold = 0.005 
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
    m[,1] = rt[i] #we might not need retention time here.
    m[,2:3] = p
    # If needed, set extra column data
    m[,1] = round(m[,1],arguments$precision_rt)
    m[,2] = round(m[,2],arguments$precision_mz)
    m[,3] = round(m[,3],0)

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
    
    # combine this function into 1 mclapply getSpectra with peaks
    MS = lapply(seq_along(CEs), function(i) getSpectra(arguments, AllIons, PEAKS, data, CEs, rt, i))
    return(list(MS, rt))
}

generateMS1forNegID <- function(arguments, MS1s, RTs, df_NegID, NegID_row){
    # Generate the MS1 for a row of NegID
    mzZoom = arguments$mzZoomWindow
    tmz = df_NegID[[1]][NegID_row,1]
    trt = df_NegID[[1]][NegID_row,2]
    MS1index = which.min(abs(RTs - trt))
    df = MS1s[[MS1index]]
    MS1_Zoom = df[ tmz - mzZoom < df[,2]
                  & df[,2] < tmz + mzZoom,, drop = FALSE]

    print( sprintf("%d, %d", NegID_row, nrow(MS1_Zoom)) )
    if(nrow(MS1_Zoom) > 0){
        mz_index = which.min(abs(MS1_Zoom[,2] - tmz))
        # this is just the closest, it's not necessarily close
        # maybe set a flag with mztol?
        if(mz_index < nrow(MS1_Zoom) ){
            dmz_scan = MS1_Zoom[(mz_index+1):nrow(MS1_Zoom),2] - MS1_Zoom[mz_index,2]
        }
    }

    MS1_Zoom = cbind(MS1_Zoom, "")
    MS1_Zoom = cbind(MS1_Zoom, "")
    MS1_Zoom = cbind(MS1_Zoom, "")
    MS1_Zoom = cbind(MS1_Zoom, df_NegID[[1]][NegID_row,3])
    MS1_Zoom = cbind(MS1_Zoom, arguments$fn_mzxml)
    MS1_Zoom = cbind(MS1_Zoom, 1)

    if(nrow(MS1_Zoom) > 0){
        MS1_Zoom[mz_index,4] = "M"
        f = df_NegID[[2]][[NegID_row]]
        MS1_Zoom[mz_index,5] = f

        if(mz_index < nrow(MS1_Zoom)){
            MS1_Zoom[(mz_index+1):(mz_index+length(dmz_scan)), 4] = get_iso_strings(arguments, dmz_scan)
        }
    }
    return( MS1_Zoom )
}

get_iso_strings <- function(arguments, dmz_scan){
    # Generate the iso string for each of the rows.
    IsoStrings = c()
    for (dmz_i in 1:length(dmz_scan)){
        isoProx = abs(arguments$Iso[[1]] - dmz_scan[dmz_i])
        isoMatches = which(isoProx < arguments$IsoThreshold)
        s = sort( isoProx[isoMatches], index.return = TRUE )
        isoMatchesSorted = isoMatches[ s$ix ]
        IsoStrings = append(IsoStrings, iso_string_genlist(arguments, isoProx, isoMatchesSorted) )
    }
    return( IsoStrings )
}

iso_string_genlist <- function(arguments, isoProx, isoMatches){
    # Generate the iso string for a single row.
    # Example output
    # 81Br(97.51%;0.00174Da);37Cl(32.40%;0.00264Da);34S(4.43%;0.00389Da);18O(0.20%;0.00455Da)
    IsoStrings = c()
    for (match_i in isoMatches){
        IsoStrings = append( IsoStrings, iso_string_gen(arguments$Iso, isoProx, match_i) )
    }
    IsoString = paste(IsoStrings, collapse=";")
    return(IsoString)
}

iso_string_gen <- function(isotopes, isoProx, match_i){
    # For a particular index, compose the isotopic string.
    return( sprintf("%s(%.2f%%;%.5fDa)", isotopes[[2]][[match_i]], isotopes[[3]][[match_i]], isoProx[[match_i]]) )
}

dfIsotopes <- function(){
    # Declare some relative masses and intensities for particular isotopes.
    # This data came from https://www.sisweb.com/mstools/isotope.htm
    arr_iso_rmass = c(0.99939,1.9958,1.00336,2.00671,1.99705,3.9941,5.99115,7.9882,2.00424,1.99795,3.99591,0.99704)
    arr_iso_symbol = c("33S","34S","13C","13C2","37Cl","37Cl2","37Cl3","37Cl4","18O","81Br","81Br2","15N")
    arr_iso_rint = c(0.7893,4.4306,1.0816,0.0117,32.3977,10.4961,3.4005,0.8501,0.2005,97.5114,48.7557,0.3613)

    return( list(arr_iso_rmass, arr_iso_symbol, arr_iso_rint) )
}

printMS1s <- function(arguments, MS1_negID, df_NegID, NegID_row){
    # Print the data to the MS1 output file.
    colnames(MS1_negID) = c("RT", "mz", "Intensity", "Isotope","Formula","PredAbundance", "Feature", "File", "Zoom")
    write.table(MS1_negID[,c("Feature", "mz", "Intensity", "Isotope","Formula","PredAbundance", "File", "Zoom"),drop=FALSE]
        , file = arguments$fn_MS1_output_wpath, row.names = FALSE, quote = FALSE
        , col.names = FALSE, append = TRUE, sep = ",")
}

calculateMS1s <- function(arguments, MS1s, RTs, df_NegID){
    # Generate the negID MS1s
    # Write MS1 data to file
    for(NegID_row in 1:nrow(df_NegID[[1]])){
        MS1_negID = generateMS1forNegID(arguments, MS1s, RTs, df_NegID, NegID_row)
        printMS1s(arguments, MS1_negID, df_NegID, NegID_row)
    }

}

printMS1header <- function(fn_MS1_output, fn_mzxml){
    # Write the header for the MS1 file
    cat(paste("Feature","m/z","Intensity","Isotope","Formula","PredAbundance"
        ,"File","Zoom", sep=","), file = fn_MS1_output, append = FALSE, sep = "\n")
}

extractCE2 <- function( arguments ){
    # Convert an mzxml file with multiple collision energies to ms2 with the fragments
    arguments$fn_mzxml_wpath = file.path(arguments$path_to_output_folder, arguments$fn_mzxml)
    AllIons <- openMSfile(arguments$fn_mzxml_wpath)  
    CEsize = getCEsize(AllIons,1)

    arguments$fn_NegID_wpath = file.path(arguments$path_to_output_folder, arguments$fn_NegID)
    df_NegID <- NegID(arguments$fn_NegID_wpath, arguments$NegID_Cols, arguments$NegID_SkipRows)
    df_NegID_mzrtid = df_NegID[,1:3]
    df_NegID_formula = df_NegID[,4]
    class(df_NegID_mzrtid) <- "numeric"
    df_NegID = list(df_NegID_mzrtid, df_NegID_formula)

    arguments$fn_MS1_output_wpath = file.path(arguments$path_to_output_folder, arguments$fn_MS1_output)
    printMS1header(arguments$fn_MS1_output_wpath, arguments$fn_mzxml)

    arguments$Iso = dfIsotopes()

    DATA = getAllSpectras(arguments, AllIons, 1, CEsize)
    MS1s = DATA[[1]]
    RTs = DATA[[2]]

    calculateMS1s(arguments, MS1s, RTs, df_NegID)
}

extractCE2( arguments )