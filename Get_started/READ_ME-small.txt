Info:
*****
# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: READ_ME

get started
***********

### requirements raw files ###
# do not change names of the raw files, structured name eg. 180616s01 
# filename must start with number, no spaces, no -,#,... symbols
# Importing data into R pipeline
- use samples after conditoning step
- paste all biological samples (samples, QCs) in folder named 'bio'
- paste 'bio' folder in 'Pipeline_metabolomics/Data/Input' folder for analysis
  ! MS: min presence of folder bio, containing at least 1 sample
  ? TROUBLESHOOTING: The required info ‘NUMBER_OF_SCANS’ in the configuration file is the first line in the example txt file. E.g. NUMBER_OF_SCANS <- 56

