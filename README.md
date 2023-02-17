# Instructions for integrating the code into the software framework:

1. Install FluoroMatch or LipidMatch from InnovativeOmics.com
2. Put Modular.r into two directories (replace existing code):
  * LipidMatch-4.2\Flow\LipidMatch_Distribution
  * LipidMatch-4.2\FluoroMatch_Modular
3. Place (genEIC.r, MS1Spectragen.r, Stats.R) into:
  * LipidMatch-4.2\Flow\LipidMatch_Distribution\LipidMatch_Libraries\Scripts
4. For the Modular version set
  * FLOW <- FALSE
  * csvInput <- FALSE #Alternatively you can keep this true and use csv inputs
  * ManuallyInputVariables <- FALSE
4. For the Flow version set
  * FLOW <- TRUE
5. Make sure to toggle the following parameters depending on your application, if both are FALSE it defaults to PFAS analysis:
  * Lipid <- TRUE
  * TWeen_pos <- FALSE
6. Make sure to follow all instructions for installing dependencies and package installation:
  * LipidMatch-4.2\FluoroMatch_Flow_Manual.docx
  * LipidMatch-4.2\FluoroMatch_Modular\FluoroMatch_Manual.docx
  * LipidMatch-4.2\FluoroMatch_Modular\Trouble_Shooting_Common_Issues.docx

If you have issues please contact Jeremy Koelmel at jeremykoelmel@gmail.com or open them up here on github