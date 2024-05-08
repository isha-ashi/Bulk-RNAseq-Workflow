sleuth.func <- function(metaData, skipLines = 0, sampleCol, factors, factorOfInt=NULL, refLevel, dataDir,
                        mappingDF, aggCol, filter_target_id = NULL, sleuth_test = "wt",
                        pvalCutoff, bValue_cutoff = 0, geneList = NULL, drop_dupGenes = FALSE, outPath){
  ## function written by Isha and Alexis to run Sleuth on outputs of kallisto. Assumes that
  ## R package "sleuth" is installed & loaded
  ## 03 January 2023
  ##
  ## metaData:    Path to the .csv file with sample information. Must have column containing
  ##              name/identifier of each sample); conditions/treatments given in other columns
  ## skipLines:   number of lines to be skipped while reading metaData file. 
  ## sampleCol:   Name of the column in metaData file containing sample names.
  ## factors:   name(s) of the column(s) in metaData containing the all the factors (multifactor) for DE analysis.
  ## factorOfInt:  name(s) of the column(s) in metaData containing the factors for which DE analysis output is needed.
  ## refLevel:    Reference level to use for the factor in factors
  ## dataDir:     filepath to the data directory 
  ## mappingDF:   data frame that stores gene annotations for each transcript/read;
  ##              must have column "target_id" containing the identifiers for each
  ##              transcript
  ## aggCol:      must be the name of one of the columns of "mappingDF"; aggregates 
  ##              targets if p-value aggregation is used
  ## filter_target_id: if non-NULL, sleuth_prep runs with target_id/transcript_id filter
  ## sleuth_test: specify the test_type to use for the analysis. "wt" for Wald test (default) or "lrt" for Likelihood Ratio test.
  ## pvalCutoff:  threshold for significance
  ## bValue_cutoff: threshold for log2 fold change. The default is 0.
  ## geneList:    (default NULL) a list of vectors, each vector containing the labels 
  ##              of genes to highlight in the volcano plots. The names of this list
  ##              should be of the form "xxxyyy" where xxx is the name of the factor 
  ##              and yyy is the variable level/group, e.g., "GroupsT1L"
  ## drop_dupGenes:  if TRUE, we drop duplicate genes from our diff. exp. lists (i.e.,
  ##              if multiple transcripts map to a gene, drop the entries for all but
  ##              one of the transcripts). Keep transcript with lowest qval.
  ## outPath:     directory to store outputs/results in 
  
  ## QC
  if(length(factors) > 1 && is.null(factorOfInt)){
    stop("factorOfInt cannot be null for multifactor analysis")
  }
  
  if(length(factors) == 1 && is.null(factorOfInt)){
    factorOfInt == factors
  }
  
  ## Create output directory
  dir.create(outPath, recursive=TRUE)
  
  ## Return results
  result <- list()
  
  # Add kallisto filepaths to sample info dataframe 
  sampleDF <- read.csv(metaData, skip = skipLines)
  colnames(sampleDF)[names(sampleDF) == sampleCol] <- "sample"
  #sampleDF[,factors] <- gsub('(.*)_\\w+', '\\1', sampleDF[,factors])
  
  ## Converting Factor columns to factor class
  sampleDF[factors] <- lapply(sampleDF[factors], factor)  ## as.factor() could also be used
  sampleDF <- sampleDF[ ,c("sample", factors)]               ## keep only "sample", factor columns
  
  ## Re-level factors 
  for (i in 1:length(factors)){  #loop thru factors by index
    fac <- factors[[i]]
    lev <- refLevel[[i]]
    sampleDF[,fac] <- relevel(sampleDF[,fac], ref = lev)
  }

  ## Getting abundance.h5 paths and adding path column to sampleDF
  abundancePaths <- c()
  for(a1 in dataDir){
    aPaths <- list.dirs(path = a1, recursive = FALSE)
    kallistoDirNames <- list.dirs(path = a1, full.names = FALSE, recursive = FALSE) ## gets directory names
    names(aPaths) <- kallistoDirNames
    if(any(kallistoDirNames %in% sampleDF$sample)){
      names(aPaths) <- kallistoDirNames
    }else{
       isPresent <- rowSums(sapply(sampleDF$sample, grepl, kallistoDirNames))
       names(isPresent) <- kallistoDirNames
       if(sum(isPresent) > 0){
          keepKallistoDirNames <- names(isPresent)[isPresent == 1]
          aPaths <- aPaths[names(aPaths) %in% keepKallistoDirNames]
          rowOrder <- max.col(sapply(sampleDF$sample, grepl, keepKallistoDirNames))
          names(aPaths) <- sampleDF$sample[rowOrder]
       }
     }

    abundancePaths <- c(abundancePaths, aPaths)
  }
  sampleDF$path <- abundancePaths[match(sampleDF$sample, names(abundancePaths))]
  
  result[["sampleDF"]] <- sampleDF  
                                                                         
  ## Create full model formula out of list of factors
  fmla <- as.formula(paste0("~", paste0(factors, collapse=" + ")))
  result[["design"]] <- fmla
  
  ## Create reduced formula out of list of factors
  if(!is.null(factorOfInt)){
    for(f1 in factorOfInt){
    redFmla <- as.formula(paste0("~", f1))
    result[["design"]] <- c(result[["design"]], redFmla)
    }
  }
    
  # Create formula object from full_model string
  # fmla <- as.formula(full_model)
  
  if(!is.null(filter_target_id)){
    sPrep <- sleuth_prep(sampleDF, read_bootstraps = FALSE, normalize = FALSE)
    filter_sleuth <- sPrep$filter_df$target_id
    filter_target_id <- filter_target_id[which(filter_target_id %in% filter_sleuth)]
  }
  
  # Create 'prepped' sleuth object
  so <- sleuth_prep(sampleDF, full_model = fmla, target_mapping = mappingDF, 
                      read_bootstrap_tpm = TRUE, filter_target_id = filter_target_id,
                      extra_bootstrap_summary = TRUE, num_cores = 1,
                      aggregation_column = aggCol, 
                      transformation_function = function(x) log2(x + 1))
  
  namesDF <- so$target_mapping[,c('gene_name','target_id')]
  
  ## Get TPM counts; add gene names to TPM count matrix
  Tpm <- sleuth_to_matrix(obj = so, which_df = "obs_raw", which_units = "tpm")
  TpmDF <- as.data.frame(Tpm)
  TpmDF_names <- merge(TpmDF, namesDF, by.x=0, by.y='target_id', all.x = TRUE)
  rownames(TpmDF_names) <- TpmDF_names$Row.names
  TpmDF_names <- subset(TpmDF_names, select=-Row.names)
  TpmDF_names <- TpmDF_names[, c(which(colnames(TpmDF_names)=="gene_name"), which(colnames(TpmDF_names)!="gene_name"))]
  
  write.csv(TpmDF_names, file.path(outPath, "tpm_names.csv"))
  result[["TPM"]] <- Tpm
  result[["TPM_names"]] <- TpmDF_names
  
  normTpm <- sleuth_to_matrix(obj = so, which_df = "obs_norm", which_units = "tpm")
  normTpmDF <- as.data.frame(normTpm)
  normTpmDF_names <- merge(normTpmDF, namesDF, by.x=0, by.y='target_id', all.x = TRUE)
  rownames(normTpmDF_names) <- normTpmDF_names$Row.names
  normTpmDF_names <- subset(normTpmDF_names, select=-Row.names)
  normTpmDF_names <- normTpmDF_names[, c(which(colnames(normTpmDF_names)=="gene_name"), which(colnames(normTpmDF_names)!="gene_name"))]
  
  result[["norm_TPM"]] <- normTpm
  result[["norm_TPM_names"]] <- normTpmDF_names
  
  ## Get normalized counts; add gene names to normalized count matrix
  normCounts <- sleuth_to_matrix(obj = so, which_df = "obs_norm", which_units = "est_counts")
  normDF <- as.data.frame(normCounts)
  normCounts_names <- merge(normDF, namesDF, by.x=0, by.y='target_id', all.x = TRUE)
  rownames(normCounts_names) <- normCounts_names$Row.names
  normCounts_names <- subset(normCounts_names, select=-Row.names)
  normCounts_names <- normCounts_names[, c(which(colnames(normCounts_names)=="gene_name"), which(colnames(normCounts_names)!="gene_name"))]
  
  write.csv(normCounts, file.path(outPath, "normalized_counts.csv"))
  write.csv(normCounts_names, file.path(outPath, "normalized_counts_names.csv"))
  
  result[["normal_counts"]] <- normCounts
  result[["normal_counts_names"]] <- normCounts_names

  
  if(sleuth_test == "wt"){
    summaryDF <- data.frame(Comparisons=character(), DE_Genes=numeric(), Up_Reg_Genes=numeric(), Down_Reg_Genes=numeric(),
                            DE_Transcripts=numeric(), Up_Reg_Transcripts=numeric(), Down_Reg_Transcripts=numeric())
    all_sigDETs <- data.frame(matrix(nrow = 0, ncol = 13))
    all_sigDEGs <- data.frame(matrix(nrow = 0, ncol = 7))
  } else { ## LRT test
    summaryDF <- data.frame(Comparison=character(), DE_Genes=numeric(), DE_Transcripts=numeric())
    all_sigDETs <- data.frame(matrix(nrow = 0, ncol = 14))
    all_sigDEGs <- data.frame(matrix(nrow = 0, ncol = 6))
  }
  
  ## Fit full model
  so <- sleuth_fit(so) # fit full model (includes all terms from 'fmla')
  
  ## Fit reduced model (only one factor)
  if(!is.null(factorOfInt)){
    so <- sleuth_fit(so, redFmla, "reduced")
  } else{ ## Fit intercept only model
    so <- sleuth_fit(obj = so, formula = ~1, fit_name = "reduced")
  }
  
  save(so,file=file.path(outPath,"sleuth_object.RData"))
  result[["sleuth_model_fit"]] <- so
  result[["sleuth_models"]] <- models(so)
  
  # Create list of factors/levels to test
  testsToRun <- colnames(so$fits$full$design_matrix)[2:ncol(so$fits$full$design_matrix)]
  result[["tests"]] <- testsToRun
  if(!is.null(factorOfInt)){
    testsToRun <- colnames(so$fits$reduced$design_matrix)[2:ncol(so$fits$reduced$design_matrix)]
    result[["tests"]] <- testsToRun
  }
    
  for(i2 in 1:length(testsToRun)){
    t1 <- testsToRun[i2]
    c1 <- gsub(factors, "", t1)
    if(!is.null(factorOfInt)){
      c1 <- gsub(factorOfInt, "", t1)
    }
    fr1 <- gsub(c1, "", t1)
    r1 <- refLevel[factors == fr1]
    if(sleuth_test == "wt"){
      # Run Wald test for each factor-level listed
      so <- sleuth_wt(so, which_beta = t1, which_model="full")
      
      ### Delete negative mean_obs values (for RIP-Seq analyses)
      mean_obs_pos <- so$tests$wt$full[[t1]][so$tests$wt$full[[t1]]$mean_obs > 0, ]
      so$tests$wt$full[[t1]] <- mean_obs_pos
      
      # Extract results from sleuth object
      deTrans <- sleuth_results(so, test = t1, test_type = "wt", which_model = "full", pval_aggregate = FALSE)
      deGenes <- sleuth_results(so, test = t1, test_type = "wt", which_model = "full", pval_aggregate = TRUE)
    }else{ ## LRT test
      # Run LRT(likelihood ratio test (LRT)) test if specifies for each factor-level listed
      
      so <- sleuth_lrt(so, null_model = "reduced", alt_model = "full") # Run a likelihood ratio test (LRT) between the two models
  
      ### Delete negative mean_obs values (for RIP-Seq analyses)
      mean_obs_pos <- so$tests$lrt$`reduced:full`[so$tests$lrt$`reduced:full`$mean_obs > 0, ]
      so$tests$lrt$`reduced:full` <- mean_obs_pos
      
      # Extract results from sleuth object
      deTrans <- sleuth_results(so, test = "reduced:full", test_type = "lrt", pval_aggregate = FALSE)
      deGenes <- sleuth_results(so, test = "reduced:full", test_type = "lrt", pval_aggregate = TRUE)
    } 
    
    # remove missing values (NA) in de_transcript data frame 
    deTrans2 <- deTrans[complete.cases(deTrans), ]
    deGenes2 <- deGenes[complete.cases(deGenes), ]
    
    # if specified, remove duplicate genes (genes w/ multiple transcripts mapped to them); keep transcript with lowest qval
    if(drop_dupGenes){ 
      deTrans2 <- deTrans2[order(deTrans2[,'gene_name'], deTrans2[,'qval']), ]
      deTrans2 <- deTrans2[!duplicated(deTrans2$gene_name), ]
      deTrans2 <- deTrans2[order(deTrans2[,'qval']), ]    #re-sort by qval
    }
      
    # Add deTrans2 to output/result
    dfName <- paste(c1, "_vs_", r1, "_allDET", sep="")
    result[[dfName]] <- deTrans2
    
    # Add deGenes2 to output/result
    dfName <- paste(c1, "_vs_", r1, "_allDEG", sep="")
    result[[dfName]] <- deGenes2
  
    # Transcripts/genes with *significant* diff. expression
    if(sleuth_test == "wt"){
      sigDET <- deTrans2[deTrans2$qval < pvalCutoff & abs(deTrans2$b) > bValue_cutoff, ]
      sigDETp <- deTrans2[deTrans2$pval < pvalCutoff & abs(deTrans2$b) > bValue_cutoff, ]
      
      ## Aggregating b-value per gene
      agg_bval <- aggregate(b ~ ensembl_gene, data = deTrans2[ ,c("ensembl_gene", "b")], sum)
      deGenes3 <- merge(deGenes2, agg_bval, by.x = "target_id", by.y = "ensembl_gene", all = FALSE, sort = FALSE) 
      sigDEG2 <- deGenes3[deGenes3$qval < pvalCutoff & abs(deGenes3$b) > bValue_cutoff, ]
    } else{
      sigDET <- deTrans2[deTrans2$qval < pvalCutoff, ]
      sigDETp <- deTrans2[deTrans2$pval < pvalCutoff, ]
      sigDEG2 <- deGenes2[deGenes2$qval < pvalCutoff, ]
    }
    
    # Create subCounts for sigDET
    subCountsResults <- c("subCountsGenes")
    if(nrow(sigDET) > 10){
      subCountsTrans <- normCounts[sigDET$target_id, ]
      subCountsResults <- c("subCountsTrans", subCountsResults)
    }
    
    # Create subCounts for sigDEG
    subCountsGenes <- normCounts_names[normCounts_names$gene_name %in% sigDEG2$gene_name, ]
    subCountsGenes <- as.matrix(subCountsGenes[,-1])
    subCountsGenes <- subCountsGenes[rowSums(subCountsGenes) > 0,]
  
    # Create/save PCA plot
    for(s1 in 1:length(subCountsResults)){
      subCounts <- get(subCountsResults[s1])
      pcaName <- sub("subCounts", "", subCountsResults[s1])
      ## browser()
      if(nrow(sigDET) > 10){
        inPCA <- as.data.frame(t(subCounts))
        inPCA <- inPCA[sampleDF$sample, ]   #put inPCA in same order as sampleDF
        ## inPCA[fr1] <- sampleDF[fr1]
        keepSamples <- sampleDF[sampleDF[,fr1] %in% c(r1, c1), "sample"] 
        inPCA <- inPCA[keepSamples,]    #keep only the samples corresp. to the reference & comparison levels
        #inPCA <- inPCA[,apply(inPCA, 2, var, na.rm=TRUE) != 0] #Exclude columns that have constant/zero values
        inPCA <- inPCA[ , !apply(inPCA, MARGIN = 2, function(x) max(x, na.rm = TRUE) == min(x, na.rm = TRUE))] #Exclude columns that have constant/zero values
        inPCA <- as.data.frame(cbind(inPCA, sampleDF[sampleDF[,fr1] %in% c(c1, r1), fr1])) 
        colnames(inPCA)[colnames(inPCA) == "sampleDF[sampleDF[, fr1] %in% c(c1, r1), fr1]"] <- fr1 
        # colnames(inPCA)[colnames(inPCA) == "sampleDF[, fr1]"] <- fr1 ###### THIS WAS THE OG LINE
        inPCA[ ,fr1] <- as.factor(inPCA[ ,fr1])
        
        samplePCA <- prcomp(inPCA[ , 1:ncol(inPCA)-1], center = TRUE, scale. = TRUE)
        # Add pca plot to result 
        pca_name <- paste(c1, "_vs_", r1, "_pca_", pcaName, sep="")
        ## pdf(paste0('/ix/djishnu/Common_Folder/',  pca_name, '.pdf'))
        pca_plot <- autoplot(samplePCA, data = inPCA, colour = fr1, size = 3) + theme(legend.position = "right", axis.text=element_text(size=12)) + geom_text(aes(label = rownames(inPCA)), nudge_y = 0.1)
        ## pdf(paste0('/ix/djishnu/Common_Folder/',  pca_name, '.pdf'))
        pdf(file.path(outPath, paste0(c1, "_vs_", r1, "_pca_", pcaName, ".pdf")))
        print(pca_plot)
        result[[pca_name]] <- print(pca_plot)
        dev.off()
        ## ggplot2::ggsave(file = file.path(outPath, paste0(c1, "_vs_", r1, "_pca_", pcaName, ".pdf")), plot = result[[pca_name]], device = "pdf")
        ## ggplot2::ggsave(file = file.path(outPath, paste0(c1, "_vs_", r1, "_pca_", pcaName, ".png")), plot = result[[pca_name]])
        
      }
    }
    
    signifFile <- file.path(outPath, paste("sig_DE_", c1, "_vs_", r1, ".xlsx", sep=""))
    sigDegName <- paste("sigDEG_", c1, "_vs_", r1, sep="")
    sigDetName <- paste("sigDET_", c1, "_vs_", r1, sep="")
    result[[sigDetName]] <- sigDET
    result[[sigDegName]] <- sigDEG2
    
    if(sleuth_test =="wt"){
      # Up and down regulated transcripts
      sigDET_Up <- sigDET[sigDET$b > 0, ]
      sigDET_Dn <- sigDET[sigDET$b < 0, ]
      sigDET_Up <- sigDET_Up[!(sigDET_Up$gene_name == ""), ] #remove the rows that has blank gene name in sigDET_Up and sigDET_Dn
      sigDET_Dn <- sigDET_Dn[!(sigDET_Dn$gene_name == ""), ] 
      
      # Add sigDEG2 to output/result
      dfName <- paste(c1, "_vs_", r1, "_allDEG_bvalue", sep="")
      result[[dfName]] <- deGenes3
      
      # Up and down regulated genes
      sigDEG2_Up <- sigDEG2[sigDEG2$b > 0, ]
      sigDEG2_Dn <- sigDEG2[sigDEG2$b < 0, ]
      sigDEG2_Up <- sigDEG2_Up[!(sigDEG2_Up$gene_name == ""), ] 
      sigDEG2_Dn <- sigDEG2_Dn[!(sigDEG2_Dn$gene_name == ""), ]
      
      # Add sigDET to result 
      sigDetName_Up <- paste(sigDetName, "_Up", sep="")
      sigDetName_Dn <- paste(sigDetName, "_Dn", sep="")
  
      # Add sigDEG to result 
      sigDegName_Up <- paste(sigDegName, "_Up", sep="")
      sigDegName_Dn <- paste(sigDegName, "_Dn", sep="")
  
      # Add (union of) all sigDET to all_sigDETs dataframe
      all_sigDETs <- rbind(all_sigDETs, sigDET) # used when performing union of PCA plot for "n" way (example n = 4 or 6)  comparisons
      
      # Add (union of) all sigDEG to all_sigDEGs dataframe
      all_sigDEGs <- rbind(all_sigDEGs, sigDEG2) # used when performing union of PCA plot for "n" way (example n = 4 or 6)  comparisons
      
      # Add a row to summary dataframe
      compStr <- paste(c1, "vs", r1)
      numDEG <- dim(sigDEG2)[1]
      numUpG <- dim(sigDEG2_Up)[1]
      numDnG <- dim(sigDEG2_Dn)[1]
      numDET <- dim(sigDET)[1]
      numUpT <- dim(sigDET_Up)[1]
      numDnT <- dim(sigDET_Dn)[1]
      
      newRow <- data.frame(Comparisons=compStr, DE_Genes=numDEG, Up_Reg_Genes=numUpG, Down_Reg_Genes=numDnG,
                           DE_Transcripts=numDET, Up_Reg_Transcripts=numUpT, Down_Reg_Transcripts=numDnT)
      
      summaryDF <- rbind(summaryDF, newRow)
      
      # Save sigDET file
      result[[sigDetName_Up]] <- sigDET_Up
      result[[sigDetName_Dn]] <- sigDET_Dn
      
      # Save sigDEG file
      result[[sigDegName_Up]] <- sigDEG2_Up
      result[[sigDegName_Dn]] <- sigDEG2_Dn
      
      write.xlsx(x=list(sigDEG=sigDEG2, sigDEG_Up=sigDEG2_Up, sigDEG_Dn=sigDEG2_Dn, 
                        sigDET=sigDET[,1:7], sigDET_Up=sigDET_Up[,1:7], sigDET_Dn=sigDET_Dn[,1:7]), file=paste0(signifFile))
      
      deResult <- c("deTrans2", "deGenes3")
      for(d1 in 1:length(deResult)){
        deResult2 <- get(deResult[d1])
        vplotName <- sub("2|3", "", sub("de", "", deResult[d1]))
        
        # Create/save volcano plot
        if (is.null(geneList)) {
          selectLab <- NULL
        } else if (is.null(geneList[[t1]])) {
          selectLab <- NULL
        } else {
          selectLab <- geneList[[t1]]
        }
        
        volc_plot <- EnhancedVolcano::EnhancedVolcano(deResult2, lab=deResult2$gene_name, x='b', y='qval',
                                                      pCutoff = pvalCutoff, pointSize = 2.0,  
                                                      xlab = bquote(~Log[2] ~ "fold change (b-value)"), axisLabSize = 11.0,
                                                      labCol= 'black', labFace = 'bold', labSize = 3.0,  
                                                      boxedLabels = TRUE, drawConnectors = TRUE, 
                                                      widthConnectors = 1.0, colConnectors = 'black',
                                                      legendLabSize = 8.0, legendIconSize = 3.0,
                                                      legendPosition = 'right', title = paste(c1, "vs", r1),
                                                      selectLab = selectLab)
        volc_plot + theme(axis.text.x = element_text(color="black", size=15, family = "Helvetica"), 
                          axis.text.y = element_text(color="black", size=15, family = "Helvetica"))
        ## ggplot2::ggsave(file = file.path(outPath, paste0(c1, "_vs_", r1, "_volcano_plot_", vplotName, ".pdf")), device = "pdf")
        ## ggplot2::ggsave(file = file.path(outPath, paste0(c1, "_vs_", r1, "_volcano_plot_", vplotName, ".png")))
        
        
        # Add volcano plot to result
        volcName <- paste(c1, "_vs_", r1, "_vPlot_", vplotName, sep="")
        result[[volcName]] <- volc_plot
      }
     
    } else {  ## LRT test
      # Add (union of) all sigDET to all_sigDETs dataframe
      all_sigDETs <- rbind(all_sigDETs, sigDET) # used when performing union of PCA plot for "n" way (example n = 4 or 6)  comparisons
      
      # Add (union of) all sigDEG to all_sigDEGs dataframe
      all_sigDEGs <- rbind(all_sigDEGs, sigDEG2) # used when performing union of PCA plot for "n" way (example n = 4 or 6)  comparisons
      
      # Add a row to summary dataframe
      compStr <- paste(c1, "vs", r1)
      numDEG <- dim(sigDEG2)[1]
      numDET <- dim(sigDET)[1]
      newRow <- data.frame(Comparison=compStr, DE_Genes=numDEG, DE_Transcr=numDET)
      summaryDF <- rbind(summaryDF, newRow)
      
      result[[sigDegName]] <- sigDEG2
      write.xlsx(x=list(sigDET=sigDET, sigDEG=sigDEG2), file=paste0(signifFile))
      }
  }
    
    result[["summaryDF"]] <- summaryDF
    result[["all_sigDETs"]] <- all_sigDETs
    
    return(result)
}

