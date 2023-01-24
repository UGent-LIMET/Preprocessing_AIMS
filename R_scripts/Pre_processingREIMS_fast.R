# Title: Software framework for the processing and statistical analysis of multivariate MS-data
# Owner: Laboratory of Integrative Metabolomics (LIMET), Ghent university
# Creator: Dr. Marilyn De Graeve
# Maintainer: <limet@ugent.be>
# Script: Part I: pre-processing REIMS


##########R Pipeline - Part I: pre-processing REIMS##########
print(Sys.time())
start_time <- Sys.time()
print("R pipeline - Part I: pre-processing REIMS fast - start!")
# Part I: pre-processing REIMS


#memory
if (CODE_RUN_MODE == CODE_DEVELOPMENT){
  memory.limit() #only check in windows possible
}


## list files from experiment 
path_data_in_bio <- file.path(path_data_in, 'bio')
setwd(path_data_in_bio)
filenames <- list.files(path_data_in_bio, pattern="*.txt", recursive = F) #not in sobfolders
amount_files <- length(filenames)

#filenames <- sample(filenames, 100) #10% min subset random files if too large for RAM, afterwrds with solve script using this output PP file

##make dir for all intermed results
dir.create(file.path(PATH, 'Data/Output', name_project, 'scans'))
path_data_out_scans <- file.path(path_data_out, 'scans')



## run over all files to find burns and collect compid matrix per file
print('loop1: find burns per sample')
files_df <- NULL
for(file_ in filenames){
  
  #file_ <- filenames[1] #for testing
  
  filename <- substr(file_, 1 ,(nchar(file_)-4))
  #print(filename)
  
  #read data
  setwd(path_data_in_bio)
  library(data.table) #fread laden
  
  #TEMP if fread fails again, slower but works:
  #scans <- read.table(file_, header=TRUE, sep="\n")
  
  try(scans <- fread(file=file_, header=TRUE, sep="\n"))  #enter seps, header= number of scans
  if(exists("scans") == FALSE){
    try(scans <- fread(file=file_, header=TRUE, sep="\n", encoding="UTF-8")) 
  }
  if(exists("scans") == FALSE){
    stop('Can not read rawt .txt files, check if UTF-16/UTF-8/ASCII format')
  }
  scans <- as.matrix(scans)
  
  #split per scan
  position_3_before_nr_scans <- unlist(gregexpr(pattern ='S', colnames(scans)))[2] #using find(S) in "NUMBER.OF.SCANS..56", 3 char befor nr
  number_of_scans <- as.numeric(substr(colnames(scans), (position_3_before_nr_scans+3), nchar(colnames(scans)))) 
  #print(number_of_scans)
  
  
  #checkpoint 1 - bad TIC chromatgram
  #only good tic gives rise to 56 and 35 scans respectively
  if(FILE_SOURCE == BURN){ #eg. sample_001.raw
    if(number_of_scans != NUMBER_OF_SCANS){
      report_bad_file <- paste0("CP1: file ", filename, " does not have a good TIC chromatogram, hence it will be excluded from the preprocessing analysis.")
      print(report_bad_file) 
      setwd(path_data_out)
      append_result_to_report(report_bad_file, paste(name_project,'_Report_bad_files_preprocessing.txt', sep=""))
      setwd(path_data_in_bio)
      #next
    }
  }
  if(FILE_SOURCE == SAMPLE){ #eg. sample.raw (had extra line: "PRECURSOR: 0 0")
    if(number_of_scans != NUMBER_OF_SCANS){
      report_bad_file <- paste0("CP1: file ", filename, " does not have a good TIC chromatogram, hence it will be excluded from the preprocessing analysis.")
      print(report_bad_file) 
      setwd(path_data_out)
      append_result_to_report(report_bad_file, paste(name_project,'_Report_bad_files_preprocessing.txt', sep=""))
      setwd(path_data_in_bio)
      #next
    }
  }
  if(FILE_SOURCE == BULK){ #eg. allinone96.raw
    #no prequisite (yet), with 96 burns/samples/measurements = 2913 # scans
  }
  
  #make empty df
  TIC_chromatogram <- data.frame(matrix(vector(), nrow=number_of_scans, ncol=(2),
                                        dimnames=list(c(), c("scan", "scan_TIC"))),
                                 stringsAsFactors=F)
  #save position of scan in scans
  position_start_scan <- data.frame(matrix(vector(), nrow=number_of_scans, ncol=(2),
                                        dimnames=list(c(), c("scannr", "line1"))),
                                 stringsAsFactors=F)
  
  #run over all scans
  line1 <- 1
  for(scannr in 1:number_of_scans){
    
    #scannr <- 1 #for testing  
    
    #set path to intermed output
    setwd(path_data_in_bio)
    #print(scannr)
    #print(line1)
    
    #save pos for each scan for later only susbet retained scans
    position_start_scan[scannr,1] <- scannr #add scan nr
    position_start_scan[scannr,2] <- line1 #add line1 = pos
    
    #amount of lines with info before mz/intensity values start
    if(FILE_SOURCE == BURN){ #eg. sample_001.raw
      lines_info <- 9
    }
    if(FILE_SOURCE == SAMPLE | FILE_SOURCE == BULK){ #eg. sample.raw (had extra line: "PRECURSOR: 0 0")
      lines_info <- 10
    }
    
    #elke scan == lines, so split
    scan <- scans[line1:(line1+lines_info),]
    scan <- as.data.frame(scan, stringsAsFactors = F)
    scanName <- scan[2,]
    #print(scanName)
    
    #mass values XXXXX X 4 BYTES, not excact but respornds to tic (more unique mz values = more likely burn) eg "MASS VALUES: 37366 x 4 BYTES"
    counts <- scan[9,]
    counts <- as.numeric(substr(counts, 14, (which(strsplit(counts, "")[[1]] == "x")-2)))
    #print(counts)
    
    #save TIC of each scan
    TIC_chromatogram[scannr,1] <- scannr #add scan nr
    TIC_chromatogram[scannr,2] <- counts #add tic
    
    line1 <- line1 + lines_info + 1
  }
  
  setwd(path_data_out_scans)
  if(FILE_SOURCE == BULK){ #eg. allinone96.raw
    #save tic-values to file, 1 file max, not per sample (too much time)
    setwd(path_data_out_scans)
    namedf <- paste0("Counts_MZ_values", filename, ".txt")
    write.table(TIC_chromatogram, file=namedf, sep ="\t", row.names = FALSE, col.names = TRUE)
    #TIC_chromatogram <- read.table(file=namedf, header=TRUE, sep="\t")
  }
  
  #Save TIC plot => no bis plot instead
  #plot(TIC_chromatogram, type="l") #show in Rstudio or save as Rplot in terminal
  #png(paste0("Counts_MZ_values", filename, ".png"), width=7, height=5, units="in", res=150)
  #plot(TIC_chromatogram, type="b")
  #dev.off()
  
  #Save TIC plot bis (ggplot) #todo ooit nice plot for fast pp
  p <- plot_Counts_MZ_values(TIC_chromatogram) 
  png(paste0("Counts_MZ_values", filename, "_bis.png"), width=7, height=5, units="in", res=150)
  plot(p)
  dev.off()  
  
  
  #checkpoint 3 - Bad TIC chromatogram
  #1st and last scan is not baseline
  #if(AMOUNT_OF_BURNS_PER_FILE == 1){ #in case of 1 burn/sample, set to 1, so skip CP3...
  #  NOISE_THRESHOLD <- .99
  #}
  region_allowed_outer_scans <- quantile(TIC_chromatogram$scan_TIC, NOISE_THRESHOLD) #10% so baseline should be lower than this
  if(TIC_chromatogram[1,2] > region_allowed_outer_scans){ #1st scan
    report_bad_file <- paste0("CP3: file ", filename, " does not have a good TIC chromatogram, hence it will be excluded from the preprocessing analysis.")
    print(report_bad_file) 
    setwd(path_data_out)
    append_result_to_report(report_bad_file, paste(name_project,'_Report_bad_files_preprocessing.txt', sep=""))
    setwd(path_data_in_bio)
    #next
  }
  if(TIC_chromatogram[nrow(TIC_chromatogram),2] > region_allowed_outer_scans){ #last scan
    report_bad_file <- paste0("CP4: file ", filename, " does not have a good TIC chromatogram, hence it will be excluded from the preprocessing analysis.")
    print(report_bad_file) 
    setwd(path_data_out)
    append_result_to_report(report_bad_file, paste(name_project,'_Report_bad_files_preprocessing.txt', sep=""))
    setwd(path_data_in_bio)
    #next
  }
  
  
  ##set threshold to select the correct burn
  #1st filter using 90% quantile (so more scans than expected burns but will be further filered below)
  THESHOLD_TIC <- quantile(TIC_chromatogram$scan_TIC, 0.90) #95%, gives mostly correct amount of burns, todo check!!! not enough, use next:
  scannr_retained_df <- TIC_chromatogram[TIC_chromatogram$scan_TIC > THESHOLD_TIC,]
  
  #order from big tic to small of top scans
  scannr_retained_df <- scannr_retained_df[order(scannr_retained_df$scan_TIC, decreasing = TRUE),] #sort from big to small tic, top from 90% quantile calc.
  scannr_retained <- scannr_retained_df[,1]
  
  #houdt rekening met afstand tsn de burns 
  #geen opeenvolgende scan mag 2e burn worden, dan ga je nr volgende hoogste tic:
  index <- 1  #index nr voor [1] als 1st elem ipv waarde v element, deze opgeslagen in top_scan
  for(top_scan in scannr_retained){
    #top_scan <- 26
    #index <- 6
    if(index < length(scannr_retained)-1){ #test until one-before-last
      if(scannr_retained[index] - 1 == scannr_retained[(index+1)] | scannr_retained[index] + 1 == scannr_retained[(index+1)]){ #eg 63 62 108 ...
        scannr_retained <- scannr_retained[-(index+1)]
      }
    }
    if(index < length(scannr_retained)-2){ 
      if(scannr_retained[index] - 1 == scannr_retained[(index+2)] | scannr_retained[index] + 1 == scannr_retained[(index+2)]){ #2 topscans actually same: eg 63 108 62 49 ...
        scannr_retained <- scannr_retained[-(index+2)]
      }
    }
    if(index < length(scannr_retained)-3){ 
      if(scannr_retained[index] - 1 == scannr_retained[(index+3)] | scannr_retained[index] + 1 == scannr_retained[(index+3)]){ #3 topscans actually same: eg 63 108 30 62 49 ...
        scannr_retained <- scannr_retained[-(index+3)]
      }#todo for bulk replace +3 by 2nd index and double loop, for sample this ok for 5burns/sample
    }
    index <- index + 1
  }
  
  #2nd filter for correct amount of scans, eg top 3 (already sorted above)
  #retain only # scans == amount of burns as defined in config file
  if(length(scannr_retained) > AMOUNT_OF_BURNS_PER_FILE){
    scannr_retained <- (scannr_retained)[1:(AMOUNT_OF_BURNS_PER_FILE)] #filter, amount scans kept dependant on config AMOUNT_OF_BURNS_PER_FILE
  }
  
  
  #checkpoint 2 - bad TIC chromatogram: 
  #if low max intensity (eg. max intens was 8000) => scannr_retained == EMPTY => skip 
  #see report in subloop retained scans + skip loading in next loop files!
  
  
  #checkpoint 5 - bad TIC chromatogram: STILL TESTING SO IN #s
  #if large drop after max burn (right to scan nr) => signal noise instead of burn => skip
  #done for each retained scan, so see below in loop
  
  #make empty df with all m/z and I of reatined scans of sample
  scans_df <-  NULL 
  
  #run over all reatained scans
  for(scannr in scannr_retained){
    #scannr <- 22 #to test
    #print('scan retained:')
    #print(scannr)

    #report is bad TIC chromatogram - (checkpoint 2 cont.)
    if(length(scannr_retained) == 0){ #if empty skip this loop
      report_bad_file <- paste0("CP2: file ", filename, " does not have a good TIC chromatogram, hence it will be excluded from the preprocessing analysis.")
      print(report_bad_file) 
      setwd(path_data_out)
      append_result_to_report(report_bad_file, paste(name_project,'_Report_bad_files_preprocessing.txt', sep=""))
      setwd(path_data_in_bio)
      next
    }
    #report is bad TIC chromatogram - (checkpoint 5 cont.)
    #if(median(scannr_retained_df[,2]) >= TIC_chromatogram$scan_TIC[scannr_retained+1]){ #if median of scans above TIC THRESHOLD (.9) >= intestiey right neigbor scan retained, to steep drop so exclude
    #  report_bad_file <- paste0("CP5: file ", filename, " does not have a good TIC chromatogram, hence it will be excluded from the preprocessing analysis.")
    #  print(report_bad_file) 
    #  setwd(path_data_out)
    #  append_result_to_report(report_bad_file, paste(name_project,'_Report_bad_files_preprocessing.txt', sep=""))
    #  setwd(path_data_in_bio)
    #  next
    #}
    
    #elke scan == amount lines, so split + only retain the wanted scan
    position_retained_scan <- position_start_scan[scannr,2]
    scan <- scans[position_retained_scan:(position_retained_scan+lines_info),]
    scan <- as.data.frame(scan, stringsAsFactors = F)
    scanName <- scan[2,]
    #print(scanName)
    
    scan <- unlist(strsplit(scan[,1], ",")) #comma seps
    scan <- as.matrix(scan)
    
    #rm first 8-9 lines (info) + last line 'scan stop' + middle line 'intensity values...'
    scan <- scan[lines_info:nrow(scan),] #rm firt info lines
    scan <- as.data.frame(scan, stringsAsFactors = F)
    scan <- scan[-nrow(scan),] #rm last line 'scan stop'
    scan <- as.data.frame(scan, stringsAsFactors = F)
    scan <- scan[-(ceiling(nrow(scan)/2)),1] #remove mid line 'intensity values'
    scan <- as.numeric(scan)
    scan <- as.data.frame(scan, stringsAsFactors = F)
    
    #now m/z and inten under each other, so split to make 2 cols
    half <- nrow(scan)/2
    #todo ev stop if not integer (instead of float)
    
    #make empty df with correct amout of rows (m/z values) and cols (samples)
    df <- data.frame(matrix(vector(), nrow=half, ncol=(2),
                            dimnames=list(c(), c(paste0(filename, "_scan", formatC(scannr, width = 3, format = "d", flag = "0"), "_MZ"), 
                                                 paste0(filename, "_scan", formatC(scannr, width = 3, format = "d", flag = "0"), "_I")))),
                     stringsAsFactors=F)
    
    #add MS values in 1st col
    df[ ,1] <- scan[1:half,1]
    
    #add intensities in 2nd col 
    df[ ,2] <- scan[(half+1):(half*2),1]
    

    #interpolate m/z
    #eg seq(50, 1200, by=0.001) as default, no binning
    if(BINNING == NO_BINNING){
      binsize <- 0.001
    }
    if(BINNING == SET_BINNING){
      binsize <- BIN_SIZE
    }
    desired_mz <- seq(START_MZ_RANGE, STOP_MZ_RANGE, by=binsize) #binsize)
    interpolated_df <- as.data.frame(approx(x=df[,1], y=df[,2], xout=desired_mz))
    colnames(interpolated_df) <- colnames(df)
    #namefile <- paste0("interpol_", name_file_scan, ".txt")
    #write.table(interpolated_df, file=namefile, sep ="\t", row.names = FALSE, col.names = TRUE)
    compid_reatined_scan <- interpolated_df
    
    
    #merge scans from 1 file
    colnames(compid_reatined_scan)[1] <- "MZ" #so merging based on non-unique m/z values
    scans_df <- merge(scans_df, compid_reatined_scan, all = TRUE) #Merging tables of different length by common columns, unqique values are NA
    
    #error first merge solve:
    if(scannr == scannr_retained[1]){
      #repeat 2nd is good
      scans_df <- merge(scans_df, compid_reatined_scan, all = TRUE) #Merging tables of different length by common columns, unqique values are NA
    }
    
    ##write matrix
    setwd(path_data_out_scans)
    #todo use write function  
    namefile <- paste0("matrix2_CopmIDs_", filename, ".txt")
    write.table(scans_df, file=namefile, sep ="\t", row.names = FALSE, col.names = TRUE)
    
    gc()
  }
} 



## make empty df for concat all samples togheter
matrix_CompIDs_all <- NULL

## run over all files 2nd time to merge indiv compid matrices into 1 matrix
print('loop2: merge samples + peak picking')
for(file_ in filenames){
  
  filename <- substr(file_, 1 ,(nchar(file_)-4))
  #print(filename)
  
  #load matrix of file
  setwd(path_data_out_scans)
  namedf <- paste0("matrix2_CopmIDs_", filename, ".txt")
  try({
    #skip if empty because of BAD chromatogram
    matrix_CompIDs <- read.table(file=namedf, header=TRUE, sep="\t")
    
    #if option AVERAGE_BURNS_PER_FILE slected, perform in M2: 
    if(MERGE_BURNS == AVERAGE_BURNS_PER_FILE){
      #average name of sample wo scannr_I
      #take average (not from mz column, start at col nr 2)
      Xfilename <- paste0('X', filename)
      matrix_CompIDs[,Xfilename] <- rowMeans(subset(matrix_CompIDs[,2:ncol(matrix_CompIDs)]), na.rm = TRUE)
      matrix_CompIDs <- matrix_CompIDs[,-c(2:(ncol(matrix_CompIDs)-1))]  #with MZ col + mean col, rm indiv column
    }
    
    matrix_CompIDs_all <- merge(matrix_CompIDs_all, matrix_CompIDs, all = TRUE) #Merging tables of different length by common columns, unqique values are NA
    
    #error first merge solve:
    if(file_ == filenames[1]){
      #repeat 2nd is good
      matrix_CompIDs_all <- merge(matrix_CompIDs_all, matrix_CompIDs, all = TRUE) 
    }
  })
  
  #head(matrix_CompIDs_all)
}  


## manual M3 merge if issues...
#setwd(path_data_out_scans)
#namedf <- "matrix3_CopmIDs_all172SCANS.txt"
#matrix_all1 <- read.table(file=namedf, header=TRUE, sep="\t")
#namedf <- "matrix3_CopmIDs_all114SCANS.txt"
#matrix_all2 <- read.table(file=namedf, header=TRUE, sep="\t")
#matrix_CompIDs_all <- cbind(matrix_all1, matrix_all2[,2:ncol(matrix_all2)]) #delete mz col

## write matrix
setwd(path_data_out_scans)
#todo use function
#use this file to for prediction!!! contains ALL, no PP/T so perfect for this; so ev todo write to input for part prediction
namedf <- paste0("matrix3_CopmIDs_all.txt")
write.table(matrix_CompIDs_all, file=namedf, sep ="\t", row.names = FALSE, col.names = TRUE)


## peak picking
#add column with sum intensities of samples
Tmatrix_CompIDs_all <- t(matrix_CompIDs_all[,2:ncol(matrix_CompIDs_all)])
sum_intensity_over_samples <-apply(Tmatrix_CompIDs_all, 2, function(x) sum(x))
matrix_CompIDs_all$sum_intensity_samples <- sum_intensity_over_samples

#meak picking using local maxima
library(data.table)
#shift lags or leads a vector by a certain amount defined as the second argument the default is to lag a vector.
#The rationale behind the below code is that each local minimum's adjucent values will be greater than itself. 
#The opposite is true for a local maximum. 
maximums <- function(x) which(x - shift(x, 1) > 0  & x - shift(x, 1, type='lead') > 0)
#minimums <- function(x) which(x - shift(x, 1) < 0  & x - shift(x, 1, type='lead') < 0)

#gives nrow of the maxima
row_nrs_to_retain <- maximums(matrix_CompIDs_all$sum_intensity_samples) #on row of sum intensity over all samples
matrix_CompIDs_all_PP <- matrix_CompIDs_all[rownames(matrix_CompIDs_all) %in% row_nrs_to_retain, ]

#write copy
namefile <- paste0("matrix3_CompIDs_all_PP.txt")
write.table(matrix_CompIDs_all_PP, file=namefile, sep ="\t", row.names = FALSE, col.names = TRUE)



## add noise threshold 
#to reduce compIDs to real compounds and not noise peaks
#threshold chosen at sum(intensities over samples) > 50 000 for same result as progenisis TODO always?
if(NOISE_REMOVAL == SET_THRESHOLD_NOISE){
  THRESHOLD_NOISE <- THRESHOLD_NOISE
}
if(NOISE_REMOVAL == DEFAULT_THRESHOLD_NOISE){
  THRESHOLD_NOISE <- quantile(matrix_CompIDs_all_PP$sum_intensity_samples, 0.95) #95%, gives mostly correct in comparison to progenmatrix, todo check!!!
  THRESHOLD_NOISE <- round(THRESHOLD_NOISE)
}

matrix_CompIDs_all_threshold <- matrix_CompIDs_all_PP[matrix_CompIDs_all_PP$sum_intensity_samples > THRESHOLD_NOISE,] 

#remove last column with sum intensity (before writing into VM)
matrix_CompIDs_all_threshold <- matrix_CompIDs_all_threshold[,-ncol(matrix_CompIDs_all_threshold)]

#write copy
namedf <- paste0("matrix3_CopmIDs_all_PP_threshold", THRESHOLD_NOISE, ".txt")
write.table(matrix_CompIDs_all_threshold, file=namedf, sep ="\t", row.names = FALSE, col.names = TRUE)



## write VM
#compID nr, add 1 rt from input and write VM for pipeline
number_compIDs <- nrow(matrix_CompIDs_all_threshold)

variableMetadata <- NULL
variableMetadata$CompID <- as.integer(seq(1, nrow(matrix_CompIDs_all_threshold), by=1))
variableMetadata$MZ <- matrix_CompIDs_all_threshold$MZ
variableMetadata$Time <- rep("1",nrow(matrix_CompIDs_all_threshold)) #put RT in mins

rest_colnames <- c("isotopes",	"adduct",	"pcgroup",	"name",	"fold",	"tstat",	"pvalue",	"mzmed",	"mzmin",	"mzmax",	"rtmed",	"rtmin",	"rtmax",	"npeaks",	"bio",	"blank"	)
rest_VM <- data.frame(matrix(vector(), nrow=number_compIDs, ncol=16,
                             dimnames=list(c(), rest_colnames)),
                      stringsAsFactors=F) #empty values for isotopes etc.
variableMetadata <- cbind(variableMetadata, rest_VM)

#add sample  intensities
variableMetadata <- cbind(variableMetadata, matrix_CompIDs_all_threshold[,2:ncol(matrix_CompIDs_all_threshold)])

#remove "X" before samplename (in P2 load added again but need to be same as VM xcms output for Part2...)
samplenames <- colnames(matrix_CompIDs_all_threshold)[2:ncol(matrix_CompIDs_all_threshold)] #retain all but first (mz column)
samplenames <- substr(samplenames, 2, nchar(samplenames))
colnames(variableMetadata) <- c("CompID", "MZ", "Time", rest_colnames, samplenames)

#save variablemetadata
setwd(path_data_in)
write_dataframe_as_txt_file(variableMetadata, paste(name_project, '_variableMetadata.txt', sep="")) #input for R_pipeline part II: statistical analysis
setwd(path_data_out)
write_dataframe_as_txt_file(variableMetadata, paste(name_project, '_variableMetadata_output_pre-processing.txt', sep="")) #copy for user



## remove all files in scan folder 
#to free up disk space, not deleted in develop mode for evaluation, ev rm as well todo
if (CODE_RUN_MODE == CODE_AUTORUN){
  #matrix1
  setwd(path_data_out_scans)
  files_rm <- list.files(path=path_data_out_scans, pattern="^matrix_")
  file.remove(files_rm) 
  
  #matrix2
  files_rm <- list.files(path=path_data_out_scans, pattern="^matrix2_")
  file.remove(files_rm) 
  
  #remark: don't delete matrix3all (use for prediction), matrix3pp, matrix3thresh and all tic png's (also if bulk tic txt)  for control
}



print("R pipeline - Part I: pre-processing REIMS fast - done!")
print(Sys.time())
end_time <- Sys.time()
print(end_time - start_time)
#
#####################
