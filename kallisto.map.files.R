kallisto.map.files <- function(inPath, rawOrTrimmed="trimmed", userName="Isha", email="idm17@pitt.edu",
                               mapTemplateFile="/ix/cigcore/utils/code/kallisto_template.sbatch",
                               transcriptsIndexPath= "/ix/cigcore/utils/mouse_indexed/mm10_kallisto.idx",
                      	       transcriptsGtfPath="/ix/cigcore/utils/mouse_anno/mm10.refGene.gtf", 
                               chromTextPath="/ix/cigcore/utils/mouse_anno/mm10.chrom.sizes", 
                               isPairedEnd=TRUE, doSubmit=TRUE, outPathMap="Auto", maxJobNum=48, hours=20){

    ## Isha wrote this function that uses kallisto_map_template.mpi to make scripts and map the reads on the cluster.
    ## INPUT:
    ## inPath: The folder that contains ORIGINAL fastq files. The cleaned fastq files are
    ## expected to be in a subfolder called "cleaned".
    ## transcriptsIndexPath: The path to the indexed transcriptome folder. Use the already indexed
    ## transcriptome.

    ## Absolute path without bash variables:
    inPath <- system(paste("echo ", inPath), intern=TRUE)
    if(unlist(strsplit(inPath,split=""))[1] !="/")
        stop("inPath must be absolute, i.e., start with '/' !")

    result <- list()
    commandS <- c()
    ## Path and files:
    splitInPath <- unlist(strsplit(inPath,split="/"))
    projPath <- paste(splitInPath[1:5], collapse = "/")
    codePath <- file.path(projPath, "code", userName)
    logPath <- file.path(codePath, "log/kallisto_align")
    dir.create(logPath, showWarnings=FALSE, recursive = TRUE)
    trimmedPath <- inPath
    if(rawOrTrimmed == "trimmed"){
      trimmedPath <- file.path(inPath,"trimmed")
    }
    print(paste("I am looking for fastq files at:", trimmedPath))
    basenames <- unlist(read.table(file = file.path(trimmedPath, "sampleList.txt")), use.names = FALSE)
    print("I found the following samples:")
    print(basenames)
    nSamples <- length(basenames)
    print(nSamples)

    ##Output path:
    if(outPathMap=="Auto"){
        outPathMap <- file.path(inPath,"mapped")
        outPathMap <- gsub(x=outPathMap, pattern="/data/", replace="/result/")
    }

    template <- readLines(mapTemplateFile)
    print(paste("isPairedEnd:", isPairedEnd))
    
    script  <- gsub(x=template, pattern = "LOGPATHOLDER", replace = logPath)
    script  <- gsub(x=script, pattern = "hoursPlaceHolder", replace = hours)
    script  <- gsub(x=script, pattern = "arrayLengthHolder", replace = nSamples-1)
    script  <- gsub(x=script, pattern = "emailHolder", replace = email)
    script  <- gsub(x=script, pattern = "dataPathPlaceHolder", replace = trimmedPath)
    script  <- gsub(x=script, pattern = "outPathPlaceHolder", replace = outPathMap)
    script <- gsub(x=script, pattern = "TRANSCRIPTSINDEXPATH", 
                   replace=transcriptsIndexPath)
    script <- gsub(x=script, pattern = "TRANSCRIPTSGTFPATH",
                   replace=transcriptsGtfPath)
    script <- gsub(x=script, pattern = "CHROMTEXTPATH",
                   replace=chromTextPath)
    
    if(!isPairedEnd){
        script<- gsub(x=script, pattern=" $DATAPATH/${FASTQNAMES[${SLURM_ARRAY_TASK_ID}]}_R2*.fastq.gz", replace="")
    } 
    
    sbatchFileName <- paste0(splitInPath[length(splitInPath)],"_kallisto_align_", rawOrTrimmed, ".sbatch")
    sbatchFilePath <- file.path(codePath, sbatchFileName)

    writeLines(script, con=sbatchFilePath)
    print(paste("I wrote:", sbatchFilePath))
    commandS <- paste("sbatch", sbatchFilePath, "\n")
    print(commandS)
    commandExit <- "exit\n"
    sleepSeconds <- 3
    if(doSubmit){
        ## Do not submit if there are too many jobs are in the queue.
        while(length(system("squeue -u $USER", intern=TRUE))-1 >= maxJobNum){
            print(paste("Still too many jobs, sleeping for",sleepSeconds/60,"minutes ...."))
            Sys.sleep(sleepSeconds)
            sleepSeconds <- sleepSeconds*1.2
        }
      myTerm <- rstudioapi::terminalCreate(show = FALSE)
      rstudioapi::terminalSend(myTerm, commandS)
      rstudioapi::terminalSend(myTerm, commandExit)
    }
    result[["commands"]] <- commandS
    mapingFile <- file.path(logPath,"mapping.RData")
    result[["mapingFile"]] <- mapingFile
    mapping <- result
    save(mapping,file=mapingFile)
    print(paste("mapping was saved in", mapingFile))
    return(mapping)
}

