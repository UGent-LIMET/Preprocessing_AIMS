## @knitr INFO
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Configuration


##########Global_settings##########

## @knitr settings
## options
RUN_CODE <- 'run this part of the pipeline'
DONT_RUN_CODE <- 'skip this part of the pipeline, keeps order: pre-processing - targeted analysis - statistical analysis - annotation'

## Adjustments
#Project name:
EXPERIMENT <- 'Preprocessing_AIMS' #structured and short, see READ_ME
POLARITY <- "negative" #{"positive", "negative"} #needed this format for annotation to work
# file_conversion and pre-processing need to be performed seperate
# only choose "both" when no pre-processing needed (eg. merge both ionisationmodes)
USER_COMMENT <- "Tutorial comment" #Add info about experiment, eg. explain (Multiple)Comparisons, to include in reports

RUN_PART_PREPROCEESSING <- RUN_CODE

#
#####################



##########Pre_processing##########
if(RUN_PART_PREPROCEESSING == RUN_CODE){
  
  ## options
  REIMS <- '.RAW folder from Waters REIMS.'
  EXACTIVE <- '.raw files from Thermo Scientific (Q)-Exactive Hybrid Quadrupole-Orbitrap Mass Spectrometer.'
  
  
  #specific options for EXACTIVE
  OPTIMIZED_PEAK_PICKING <- 'Optimisation of peak picking parameters with IPO library, very slow.' 
  FIXED_PEAK_PICKING <- 'Use default values of slides G4M software, faster and similar results as optimized.'
  # Choose FIXED_PEAK_PICKING for targeted analysis, gives bit less compounds and has wider rt-range so easier retrieve standards (Targeted_approach)
  
  SAVE_PLOTS <- 'for optimalisatiton purposes: see plots when working in development phase (@Rstudio, needs graphical interphase to work). slower.'
  DONT_SAVE_PLOTS <- 'Without saving plots in output. faster.'
  
  
  #specific options for REIMS
  BURN <- 'output after processing Progenisis Bridge, 1 raw file contains 1 cleaned up burn.'
  SAMPLE <- 'output after MassLynx, 1 raw file contains 1 sample (with 1 or more burns).'
  BULK <- 'output after MassLynx, 1 raw file contains multiple sample (with 1 or more burns each).'
  
  #MERGE_BURNS options
  INDIVIDUAL_BURNS_PER_FILE <- 'keeps intensity of each burn present per 1 raw file'
  AVERAGE_BURNS_PER_FILE <- 'takes the average intensity of burns present per 1 raw file'
  
  #BINNING options
  NO_BINNING <- "no binning, MZs are evaluated as come from aquisition"
  SET_BINNING <- "set manual bin size using BIN_SIZE"
  
  #NOISE_REMOVAL options
  DEFAULT_THRESHOLD_NOISE <- "default noise removal based on distriburtion intensities peaks (retain 95% quantile)"
  SET_THRESHOLD_NOISE <- "set manual noise removal insensity"

  
  
  ## Adjustments
  INSTRUMENT <- REIMS 
  
  
  #specific adjustments for EXACTIVE
  PARAMETERS_PREPROCESSING <- FIXED_PEAK_PICKING  #parameters for peak picking, rt align and filling
  
  PLOTTING <- SAVE_PLOTS
  
  
  #specific adjustments for REIMS
  FILE_SOURCE <- SAMPLE #SAMPLE, BURN or BULK, see above explanation
  
  AMOUNT_OF_BURNS_PER_FILE <- 1 #how may burns per raw file, needs integer eg. 1, 2, 3
  MERGE_BURNS <- INDIVIDUAL_BURNS_PER_FILE
  
  NUMBER_OF_SCANS <- 114 #how many timepoints are calculated per file, check fisrt run sample .txt file "NUMBER.OF.SCANS..56" eg. 35 (bridge output), 56 (~30s), 114 (~1min)
  
  NOISE_THRESHOLD <- 0.90 #ratio max burn intensity / baseline intensity: 0.90 default, 0.95 or higher needed in case of high baseline
  
  START_MZ_RANGE <- 50    #settings dat-aquisition eg. from 50-1200 Da: START_MZ_RANGE <- 50
  STOP_MZ_RANGE <- 1200   #settings dat-aquisition eg. from 50-1200 Da: STOP_MZ_RANGE <- 1200
  
  BIN_SIZE <- 0.1
  BINNING <- NO_BINNING
  
  THRESHOLD_NOISE <- 50000
  NOISE_REMOVAL <- DEFAULT_THRESHOLD_NOISE
  
}
#
#####################
