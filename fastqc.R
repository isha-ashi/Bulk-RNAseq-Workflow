fastqc <- function(dataPath, toExtractFastaFile=NULL, rawOrTrimmed="raw", userName="Isha",
                   sbatchTemplateFile="/ix/cigcore/utils/code/fastqc_template.sbatch",
                   doRemoveOriginalFastq=FALSE, doSubmit=TRUE, email="idm17@pitt.edu",
                   suffix="fastq.gz", charNumAfterBasename="Auto", maxJobNum=48, doRemoveBam=FALSE,
                   doWrite=TRUE, transferPrefix=NULL, minlen=NULL, hours=20){
  library(stringr)
  library(stringi)
  
  ## Isha wrote this function that uses clean_fastq.sh on the cluster to clean all raw fastq files.
  ## in a given directory. This function can be used by clean_fastqS.R.
  ## INPUT:
  ## dataPath: The path to the folder that contains fastq files.
  ## charNumAfterBasename: Count the number of characters in the file name from last including "?R1".
  ## doRrna: If TRUE, rRNA will be removed. Set to TRUE for RNA, and FALSE for DNA, data.   
  print(paste("Cleaning fastq files at", dataPath))
  
  ## QC
  if(unlist(strsplit(dataPath,split=""))[1] !="/")
    stop("dataPath must be absolute, i.e., start with '/' !")
  
  result <- list()
  commandS <- c()
  
  ## Data and Path:
  qcOutPath <- file.path(dataPath, "qc", rawOrTrimmed)
  if(rawOrTrimmed == "trimmed"){
    dataPath <- file.path(dataPath, "trimmed")
  }
  splitDataPath <- unlist(strsplit(dataPath,split="/"))
  projPath <- paste(splitDataPath[1:5], collapse = "/")
  codePath <- file.path(projPath, "code", userName)
  logPath <- file.path(codePath, "log/fastqc")
  dir.create(logPath, showWarnings=FALSE, recursive = TRUE)
  
  ## Files:
  files <- list.files(dataPath, pattern = suffix)
  
  ## Number of characters in the file names after the basename.
  if(charNumAfterBasename=="Auto"){
    if(suffix=="fastq.gz"){
      if(grepl("\\R1", files[1])|grepl("\\R2", basename(files[1]))){
        print("R1/R2 exist in file names!")
        
      } else{
        stop("File names don't include R1/R2!!")
      }
      if(lengths(regmatches(files[1], gregexpr("_R", files[1]))) >1){
        stop("File names should not have more than one R1/R2 in their names!!!")
      }
      charNumAfterBasename <-
        nchar(tail(unlist(strsplit(files[1], "\\_R")), n=1))+nchar("_R")
      } else
      if(suffix=="bam"){
        charNumAfterBasename <- 4
      } else {
        stop("charNumAfterBasename could not be determined automatically!")
      }
  }
  print(paste("charNumAfterBasename:", charNumAfterBasename))
  basenames <- unique(gsub(paste0('.{',charNumAfterBasename,'}$'), '', files))
  cat(basenames, file = file.path(dataPath, "sampleList.txt"), sep = "\n")
  print("I found the following samples:")
  print(basenames)
  nSamples <- length(basenames)
  print(nSamples)
  filEndpattern <- paste0("_R", unique(gsub(paste0('.*_R'), '', files)))
  print(filEndpattern)
  result[["basenames"]] <- basenames
  
  template <- readLines(sbatchTemplateFile)
  ## Common parameters:
  doBam2fastq <- suffix=="bam"
  if(!is.null(minlen)){
    minlenReplace <- paste("-l", minlen)
  } else {
    minlenReplace <- ""
  }
  
  fastqcJob <- paste0("fastqc -o $QCOUTPATH $DATAPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}", 
                       filEndpattern)
  
  print("Preparing a script and submitting for:")
  template <- readLines(sbatchTemplateFile)
  script  <- gsub(x=template, pattern = "LOGPATHOLDER", replace = logPath)
  script  <- gsub(x=script, pattern = "hoursPlaceHolder", replace = hours)
  script  <- gsub(x=script, pattern = "arrayLengthHolder", replace = nSamples-1)
  script  <- gsub(x=script, pattern = "emailHolder", replace = email)
  script  <- gsub(x=script, pattern = "dataPathPlaceHolder", replace = dataPath)
  script  <- gsub(x=script, pattern = "minlenPlaceHolder", replace = minlenReplace)
  script  <- gsub(x=script, pattern = "qcFolderNameHolder", replace = qcOutPath)
  script  <- gsub(x=script, pattern = "jobCommandHolder1", replace = fastqcJob[1])
  if(length(fastqcJob) == 2){
    script <- gsub(x=script, pattern = "jobCommandHolder2", replace = fastqcJob[2])
  } else{
    script <- gsub(x=script, pattern = "jobCommandHolder2", replace = "")
  }
  
  sbatchFileName <- paste0(splitDataPath[length(splitDataPath)],"_fastqc_", rawOrTrimmed, ".sbatch")
  sbatchFilePath <- file.path(codePath, sbatchFileName)
  
  if(doWrite){
    writeLines(script, con=sbatchFilePath)
    print(paste("I wrote:", sbatchFilePath))
    commandS <- paste0("sbatch ",sbatchFilePath, "\n")
    print(commandS)
    sleepSeconds <- 3
    if(doSubmit){
      ## Do not submit if there are too many jobs are in the queue.
      while(length(system("squeue -u $USER", intern=TRUE))-1 >= maxJobNum){
        print(paste("Still too many jobs (",maxJobNum,"), sleeping for",
                    sleepSeconds/60,"minutes ...."))
        Sys.sleep(sleepSeconds)
        sleepSeconds <- sleepSeconds*1.2
      }
      myTerm <- rstudioapi::terminalCreate(show = FALSE)
      rstudioapi::terminalSend(myTerm, commandS)
      Sys.sleep(1)
      repeat{
          Sys.sleep(0.1)
          if(rstudioapi::terminalBusy(myTerm) == FALSE){
             print("Code Executed")
             break
           }
      }
      rstudioapi::terminalKill(myTerm)
      print(paste("Submitted at:", Sys.time()))
    }
  }
  result[["commands"]] <- commandS

  ## How to copy QC to local?
  if(!is.null(transferPrefix)){
    toPath <- gsub(dataPath, pattern="^.*/proj/", replacement="~/OneDrive\ -\ University\ of\ Pittsburgh/proj")
    transferCommand <- paste0("scp ", transferPrefix, file.path(dataPath,"qc "), toPath)
    print("Copy QC reports: Run the following from the local computer AFTER the cleaning jobs are done.")
    print(transferCommand)
    result[["transferCommand"]] <- transferCommand
  }
  cleaning <- result
  save(cleaning,file=file.path(logPath,"cleaning.RData"))
  return(cleaning)
}