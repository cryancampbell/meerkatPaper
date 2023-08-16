## using counts and metadata, generate residuals for emmreml analysis

# test a change

# batchList <- "allBatches"
# focalTreatment <- "LPS"
# minGenes <- 10000
# minLogCPM <- 5
# 
# Kmat <- TRUE
# Kfile <- "~/Dropbox (Personal)/meerkats/KMPsamples/KMPrelatedness_matrix.RDS"
# 
# singleTrt <- FALSE
# combineTrt <- FALSE
# 
# if (singleTrt) {
#   gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"only_",minGenes,"genes_",minLogCPM,"lCPM")
# } else if (combineTrt)  {
#   gatherSubset <- paste0(batchList,"_",focalSex,"_allTrt_",minGenes,"genes_",minLogCPM,"lCPM")
# } else {
#   gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"_",minGenes,"genes_",minLogCPM,"lCPM")
# }
# 
# perms <- TRUE
# permVersion <- 1
# permVar <- "treatment"
# 
# 
# focalTreatment <- "Gard"; modelName <- "G"
# focalSex <- "F"; LM <- "basic+Preg+WeightwKmatrix"
#focalSex <- "M"; LM <- "basic+WeightwKmatrix"



runEmreml <- function(workDir, gatherSubset, model, LM, perms = FALSE, permVar = "treatment", permNum = 1, 
                      Kmat = TRUE, Kfile = "~/Dropbox (Personal)/meerkats/KMPsamples/KMPrelatedness_matrix.RDS",
                      permVersion = 1){
  dataDir <- paste0(workDir,"data/")
  outDir <- paste0(workDir,"output/")
  figDir <- paste0(workDir,"figures/")
  
  dir.create(paste0(figDir,gatherSubset,"/model",model))
  modelDir <- paste0(figDir,gatherSubset,"/model",model)
  if (perms == TRUE) {
    modelDir <- paste0(figDir,gatherSubset,"/model",model,"/perms")
    dir.create(modelDir)
  }
  
  #read in data
  subsetMeta <- readRDS(file = paste0(outDir,"/",gatherSubset,"_metadata.RDS"))
  subsetCountsDF <- readRDS(file = paste0(outDir,"/",gatherSubset,"_counts.RDS"))
  filterCountsDF <- readRDS(file = paste0(outDir,"/",gatherSubset,"_filteredCounts.RDS"))
  
  resids <- readRDS(file = paste0(outDir,"/",gatherSubset,"_model",model,"_resids.RDS"))
  
  if (Kmat) {
    K <- readRDS(Kfile)
  }
  
  ### load counts, residuals, and metadata
  #subsetMeta <- read.table(file = paste0(outDir,modelName,"/",modelName,"_metadata"))
  #filterCountsDF <- read.table(file = paste0(outDir,modelName,"/",modelName,"_filtCounts"))
  colnames(filterCountsDF) <- subsetMeta$sampleID
  
  ### CHECK THAT THEY MATCH?!?
  
  #resids <- read.table(file = paste0(outDir,modelName,"/",modelName,"_residuals"))
  resids <- as.matrix(resids)
  colnames(resids) <- subsetMeta$sampleID
  
  nSamples <- dim(subsetMeta)[1]
  nIndividuals <- length(unique(subsetMeta$correctsampleID))
  nGenes <- dim(resids)[1]
  
  
  #Z - individual to sample mapping
  sampleIndMap <- diag(0,nrow=nSamples,ncol=nIndividuals)
  rownames(sampleIndMap) <- subsetMeta$sampleID
  colnames(sampleIndMap) <- unique(subsetMeta$correctsampleID)
  
  for (i in 1:nSamples) {
    for (j in 1:nIndividuals) {
      sampleIndMap[i,j] <- 1 * (subsetMeta$correctsampleID[i] == colnames(sampleIndMap)[j])
    }
  }
  
  if (Kmat) {
    #start with identity
    identityK <- diag(1,nrow=nIndividuals,ncol=nIndividuals)
    rownames(identityK) <- unique(subsetMeta$correctsampleID)
    colnames(identityK) <- unique(subsetMeta$correctsampleID)
    
    Ksubset <- K[rownames(K) %in% subsetMeta$kmpID,colnames(K) %in% subsetMeta$kmpID]
    
    #look through identityK and replace 0's with correct value from Ksubset
    for (i in 1:nIndividuals) {
      for (j in 1:nIndividuals) {
        if (i != j) {
          iKMPid <- subset(subsetMeta, correctsampleID == rownames(identityK)[i])$kmpID[1]
          jKMPid <- subset(subsetMeta, correctsampleID == rownames(identityK)[j])$kmpID[1]
          
          if (iKMPid %in% rownames(Ksubset) & jKMPid %in% colnames(Ksubset)) {
            identityK[i,j] <- Ksubset[which(rownames(Ksubset) == iKMPid),
                                      which(colnames(Ksubset) == jKMPid)]
          }
        }
      }
    }
    
  } else {
    #K relatedness, none as of yet
    identityK <- diag(1,nrow=nIndividuals,ncol=nIndividuals)
  }
  
  
  subsetMeta$permID <- paste0(subsetMeta$animal_id,"_",subsetMeta$captureref)
  
  fixed_cov <- subsetMeta[which(colnames(subsetMeta) %in% c("Trt","sex","status_dom_sub","ageAtSample","permID","animal_id","pregLact","climSeason","avgWeight6mo","captureref"))]
  
  colnames(fixed_cov)[which(colnames(fixed_cov) == "Trt")] <- "treatment"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "status_dom_sub")] <- "status"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "ageAtSample")] <- "age"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "pregLact")] <- "reprod"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "climSeason")] <- "season"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "avgWeight6mo")] <- "weight"
  
  #singleTrt <- FALSE
  #combineTrt <- FALSE
  ## learn if this is single treatment from the data
  if (length(unique(fixed_cov$treatment)) == 1) {
    singleTrt <- TRUE
    trtName <- unique(fixed_cov$treatment)
  } else if (length(unique(fixed_cov$treatment)) == 2) {
    trtName <- unique(fixed_cov$treatment)[which(unique(fixed_cov$treatment) != "Ctrl")]
  } else {
    combineTrt <- TRUE
    trtName <- "all"
  }

  fixed_cov$age <- scale(fixed_cov$age, scale = F, center = T)
  fixed_cov$weightUnScaled <- fixed_cov$weight
  animalWeights <- unique(fixed_cov$weight)
  scaledWeights <- scale(animalWeights, scale = F, center = T)
  ###WORKTODO###
  #need to fix the weights to be truly scaled (each individual counts once)
  for (n in 1:dim(fixed_cov)[1]) {
    matchRow <- which(animalWeights == fixed_cov$weightUnScaled[n])
    fixed_cov$weight[n] <- scaledWeights[matchRow]
  }
  
  #same thing for status only
  dWeights <- unique(subset(fixed_cov, status == "D")$weightUnScaled)
  dWeightsScaled <- scale(dWeights, scale = F, center = T)
  sWeights <- unique(subset(fixed_cov, status == "S")$weightUnScaled)
  sWeightsScaled <- scale(sWeights, scale = F, center = T)
  fixed_cov$weightStatusScaled <- 0
  
  for (n in 1:dim(fixed_cov)[1]) {
    if (fixed_cov$status[n] == "D") {
      matchRow <- which(dWeights == fixed_cov$weightUnScaled[n])
      fixed_cov$weightStatusScaled[n] <- dWeightsScaled[matchRow]
    } else {
      matchRow <- which(sWeights == fixed_cov$weightUnScaled[n])
      fixed_cov$weightStatusScaled[n] <- sWeightsScaled[matchRow]
    }
  }
  
  if (singleTrt == FALSE & combineTrt == FALSE) {
    fixed_cov$treatment <- factor(fixed_cov$treatment,
                                  levels = c("Ctrl",
                                           trtName))
  } else if (singleTrt == FALSE & combineTrt == TRUE) {
    fixed_cov$treatment <- factor(fixed_cov$treatment,
                                  levels = c("Ctrl",
                                             "Dex",
                                             "Gard",
                                             "LPS"))
  }
  
  fixed_cov$status <- factor(fixed_cov$status, levels = c("S","D"))
  fixed_cov$sex <- factor(fixed_cov$sex, levels = c("F","M"))
  fixed_cov$reprod <- factor(fixed_cov$reprod, levels = c("none","preg","lact"))
  fixed_cov$preg <- ifelse(fixed_cov$reprod == "preg", TRUE, FALSE)
  
  fixed_cov$season <- factor(fixed_cov$season, levels = c("warm-wet","cold-dry"))
  
  if ( perms ) {
    if (permVersion == 1) {
      trueFixed_Cov <- fixed_cov
      if ( permVar == "treatment" ) {
        #permute treatment variable here
        possibleTrts <- unique(fixed_cov$treatment)
        
        chall <- unique(fixed_cov$permID)[22]
        for (chall in unique(fixed_cov$permID)) {
          covRows <- which(fixed_cov$permID == chall)
          
          actualTrts <- fixed_cov$treatment[covRows]
          #permTrts <- actualTrts
          
          #if there are 2 entries, randomly flip or not
          if (length(actualTrts) == 2) {
            if (runif(1) < .5) {
              permTrts <- actualTrts[c(2,1)]
            } else {
              permTrts <- actualTrts
            }
          } else {
            #otherwise randomly pick trt or control
            if (runif(1) < .5) {
              permTrts <- possibleTrts[1]
            } else {
              permTrts <- possibleTrts[2]
            }
          }
          fixed_cov$treatment[covRows] <- permTrts
        }
        
      } else if ( permVar %in% c("status","preg","weight","age") ) {
        #fixed_cov$animalStatus <- paste0(fixed_cov$animal_id,"_",fixed_cov$status)
        #get all captures and their status/weight/age/preg
        captStat <- unique(fixed_cov[, which(colnames(fixed_cov) %in% c("captureref",permVar))])
        varColNum <- which(colnames(fixed_cov) == permVar)
        trueCaptStat <- captStat
        
        #permute statuses
        captStat[,2] <- sample(captStat[,2], size = length(captStat[,2]), replace = FALSE)
        
        #fill new status back into fixed_cov DF
        for (sampleEntry in 1:dim(fixed_cov)[1]) {
          captID <- fixed_cov$captureref[sampleEntry]
          permStatus <- subset(captStat, captureref == captID)[,2]
          fixed_cov[sampleEntry,varColNum] <- permStatus
        }
      } 
    } else if (permVersion == 2) {
      #shuffle order of columns in residual object
      trueResids <- resids
      resids <- resids[,sample(1:dim(resids)[2], size = dim(resids)[2], replace = FALSE)]
      
      ### shuffle WITH BLOCKS?!?
      
    } else if (permVersion == 3) {
      trueFixed_Cov <- fixed_cov
      #get all captures and their status+weight+age+preg
      #### TO-DO sort out male/female here (drop preg in male models)
      captStat <- unique(fixed_cov[, which(colnames(fixed_cov) %in% c("captureref","status","weight","age","preg"))])
      #varColNum <- which(colnames(fixed_cov) == permVar)
      trueCaptStat <- captStat
      
      #permute CAPTURES
      captStat[,1] <- sample(captStat[,1], size = length(captStat[,1]), replace = FALSE)
      
      #fill new status/weight/age/preg back into fixed_cov DF
      for (sampleEntry in 1:dim(fixed_cov)[1]) {
        captID <- fixed_cov$captureref[sampleEntry]
        
        fixed_cov$status[sampleEntry] <- subset(captStat, captureref == captID)$status
        fixed_cov$age[sampleEntry] <- (subset(captStat, captureref == captID)$age)
        fixed_cov$weight[sampleEntry] <- (subset(captStat, captureref == captID)$weight)
        fixed_cov$preg[sampleEntry] <- subset(captStat, captureref == captID)$preg
      }
    }
  }
  
  #?#?# differentiate trt-ctrl v condition analysis here
  if (singleTrt) {
    if (LM == "nested") {
      LM <- "basic"
    }
  }
  
  
  if ( LM == "basic" ) {
    
    if (singleTrt) {
      
      #basic model
      fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                   fixed_cov$age)
      
      emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
      rownames(emmremlOutput) <- rownames(resids)
      colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age',
                                   'var_beta_intercept','var_beta_status','var_beta_age',
                                   'pval_intercept','pval_status','pval_age')
    } else {

      #basic model
      fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                   fixed_cov$status + 
                                   fixed_cov$age)
      
      emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
      rownames(emmremlOutput) <- rownames(resids)
      colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_status','beta_age',
                                   'var_beta_intercept','var_beta_treatment','var_beta_status','var_beta_age',
                                   'pval_intercept','pval_treatment','pval_status','pval_age')
    }
      } else if ( LM == "basic+Weight" ) {
        
        if (singleTrt) {
          
          #basic model
          fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                       fixed_cov$age +
                                       fixed_cov$weight)
          
          emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
          rownames(emmremlOutput) <- rownames(resids)
          colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_weight',
                                       'var_beta_intercept','var_beta_status','var_beta_age','var_beta_weight',
                                       'pval_intercept','pval_status','pval_age','pval_weight')
        } else {
          
          #basic model
          fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                       fixed_cov$status + 
                                       fixed_cov$age +
                                       fixed_cov$weight)
          
          emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
          rownames(emmremlOutput) <- rownames(resids)
          colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_status','beta_age','beta_weight',
                                       'var_beta_intercept','var_beta_treatment','var_beta_status','var_beta_age','var_beta_weight',
                                       'pval_intercept','pval_treatment','pval_status','pval_age','pval_weight')
        }
        
        } else if ( LM == "basicMF+Weight" ) {
          
          if (singleTrt) {
            
            #basic model
            fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                         fixed_cov$age +
                                         fixed_cov$weight +
                                         fixed_cov$sex)
            
            emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
            rownames(emmremlOutput) <- rownames(resids)
            colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_weight','beta_sex',
                                         'var_beta_intercept','var_beta_status','var_beta_age','var_beta_weight','var_beta_sex',
                                         'pval_intercept','pval_status','pval_age','pval_weight','pval_sex')
          } else {
            
            #basic model
            fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                         fixed_cov$status + 
                                         fixed_cov$age +
                                         fixed_cov$weight +
                                         fixed_cov$sex)
            
            emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
            rownames(emmremlOutput) <- rownames(resids)
            colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_status','beta_age','beta_weight','beta_sex',
                                         'var_beta_intercept','var_beta_treatment','var_beta_status','var_beta_age','var_beta_weight','var_beta_sex',
                                         'pval_intercept','pval_treatment','pval_status','pval_age','pval_weight','pval_sex')
          }
          
        }else if (LM == "nested") {
        
        fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                     fixed_cov$status:fixed_cov$treatment + 
                                     fixed_cov$age)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_age',
                                     'beta_Ctrl:Status','beta_Trt:Status',
                                     'var_beta_intercept','var_beta_treatment','var_beta_age',
                                     'var_beta_Ctrl:Status','var_beta_Trt:Status',
                                     'pval_intercept','pval_treatment','pval_age',
                                     'pval_Ctrl:Status','pval_Trt:Status')
        
        #want to compare S - D in each of LPS and Ctrl
        fixedCovMM[,4] <- ifelse(fixed_cov$treatment == trtName, 0,
                                 ifelse(fixed_cov$status == "D", 1, -1))
        
        fixedCovMM[,5] <- ifelse(fixed_cov$treatment == "Ctrl", 0,
                                 ifelse(fixed_cov$status == "D", 1, -1))
      } else if (LM == "nested+Weight") {
        
        fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                     fixed_cov$status:fixed_cov$treatment + 
                                     fixed_cov$age +
                                     fixed_cov$weight)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_age','beta_weight',
                                     'beta_Ctrl:Status','beta_Trt:Status',
                                     'var_beta_intercept','var_beta_treatment','var_beta_age','var_beta_weight',
                                     'var_beta_Ctrl:Status','var_beta_Trt:Status',
                                     'pval_intercept','pval_treatment','pval_age','pval_weight',
                                     'pval_Ctrl:Status','pval_Trt:Status')
        
        #want to compare S - D in each of LPS and Ctrl
        fixedCovMM[,5] <- ifelse(fixed_cov$treatment == trtName, 0,
                                 ifelse(fixed_cov$status == "D", 1, -1))
        
        fixedCovMM[,6] <- ifelse(fixed_cov$treatment == "Ctrl", 0,
                                 ifelse(fixed_cov$status == "D", 1, -1))
      } else if (LM == "nestedWPreg") {
        
        fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                     fixed_cov$status:fixed_cov$treatment + 
                                     fixed_cov$age + 
                                     fixed_cov$preg)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_age','beta_preg',
                                     'beta_Ctrl:Status','beta_Trt:Status',
                                     'var_beta_intercept','var_beta_treatment','var_beta_age','var_beta_preg',
                                     'var_beta_Ctrl:Status','var_beta_Trt:Status',
                                     'pval_intercept','pval_treatment','pval_age','pval_preg',
                                     'pval_Ctrl:Status','pval_Trt:Status')
        
        #want to compare S - D in each of LPS and Ctrl
        fixedCovMM[,5] <- ifelse(fixed_cov$treatment == trtName, 0,
                                 ifelse(fixed_cov$status == "D", 1, -1))
        
        fixedCovMM[,6] <- ifelse(fixed_cov$treatment == "Ctrl", 0,
                                 ifelse(fixed_cov$status == "D", 1, -1))
      } else if (LM == "nested+Preg+Weight") {
        fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                     fixed_cov$status:fixed_cov$treatment + 
                                     fixed_cov$age + 
                                     fixed_cov$preg + 
                                     fixed_cov$weight)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_age','beta_preg','beta_weight',
                                     'beta_Ctrl:Status','beta_Trt:Status',
                                     'var_beta_intercept','var_beta_treatment','var_beta_age','var_beta_preg','var_beta_weight',
                                     'var_beta_Ctrl:Status','var_beta_Trt:Status',
                                     'pval_intercept','pval_treatment','pval_age','pval_preg','pval_weight',
                                     'pval_Ctrl:Status','pval_Trt:Status')
        
        #want to compare S - D in each of LPS and Ctrl
        fixedCovMM[,6] <- ifelse(fixed_cov$treatment == trtName, 0,
                                 ifelse(fixed_cov$status == "D", 1, -1))
        
        fixedCovMM[,7] <- ifelse(fixed_cov$treatment == "Ctrl", 0,
                                 ifelse(fixed_cov$status == "D", 1, -1))
        
      } else if (LM == "nestedPreg") {
        
        ###only test preg-v-non here
        fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                     fixed_cov$preg:fixed_cov$status + 
                                     fixed_cov$age)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age',
                                     'beta_Sub:Preg','beta_Dom:Preg',
                                     'var_beta_intercept','var_beta_status','var_beta_age',
                                     'var_beta_Sub:Preg','var_beta_Dom:Preg',
                                     'pval_intercept','pval_status','pval_age',
                                     'pval_Sub:Preg','pval_Dom:Preg')
        
        #need to center the nested variables, -1,0,1
        #dom + preg
        fixedCovMM[,4] <- ifelse(fixed_cov$status == "S" & fixed_cov$preg == TRUE, 1, 
                                 ifelse(fixed_cov$status == "S" & fixed_cov$preg == FALSE, -1, 0))
        #sub + preg
        fixedCovMM[,5] <- ifelse(fixed_cov$status == "D" & fixed_cov$preg == TRUE, 1, 
                                 ifelse(fixed_cov$status == "D" & fixed_cov$preg == FALSE, -1, 0))
        
      } else if (LM == "nestedStatus") {
        
        ###only test preg-v-non here
        fixedCovMM <- model.matrix(~ fixed_cov$preg + 
                                     fixed_cov$status:fixed_cov$preg + 
                                     fixed_cov$age)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_preg','beta_age',
                                     'beta_NonPreg:Status','beta_Preg:Status',
                                     'var_beta_intercept','var_beta_preg','var_beta_age',
                                     'var_beta_NonPreg:Status','var_beta_Preg:Status',
                                     'pval_intercept','pval_preg','pval_age',
                                     'pval_NonPreg:Status','pval_Preg:Status')
        
        #need to center the nested variables, -1,0,1
        #non-preg + dom -> -1
        fixedCovMM[,4] <- ifelse(fixed_cov$preg == TRUE, 0, 
                                 ifelse(fixed_cov$status == "D", 1, -1))
        #preg + dom -> -1
        fixedCovMM[,5] <- ifelse(fixed_cov$preg == FALSE, 0, 
                                 ifelse(fixed_cov$status == "D", 1, -1))
        
      } else if (LM == "basic+Preg") {
        if (singleTrt) {
          
          fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                       fixed_cov$age +
                                       fixed_cov$preg)
          
          emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
          rownames(emmremlOutput) <- rownames(resids)
          colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg',
                                       'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg',
                                       'pval_intercept','pval_status','pval_age','pval_preg')
        } else {
          fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                       fixed_cov$status + 
                                       fixed_cov$age +
                                       fixed_cov$preg)
          
          emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
          rownames(emmremlOutput) <- rownames(resids)
          colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_status','beta_age','beta_preg',
                                       'var_beta_intercept','var_beta_treatment','var_beta_status','var_beta_age','var_beta_preg',
                                       'pval_intercept','pval_treatment','pval_status','pval_age','pval_preg')
          
        }
      } else if (LM == "basic+Preg+Weight") {
        if (singleTrt) {
          
          fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                       fixed_cov$age +
                                       fixed_cov$preg +
                                       fixed_cov$weight)
          
          emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
          rownames(emmremlOutput) <- rownames(resids)
          colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg','beta_weight',
                                       'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg','var_beta_weight',
                                       'pval_intercept','pval_status','pval_age','pval_preg','pval_weight')
        } else {
          fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                       fixed_cov$status + 
                                       fixed_cov$age +
                                       fixed_cov$preg +
                                       fixed_cov$weight)
          
          emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
          rownames(emmremlOutput) <- rownames(resids)
          colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_status','beta_age','beta_preg','beta_weight',
                                       'var_beta_intercept','var_beta_treatment','var_beta_status','var_beta_age','var_beta_preg','var_beta_weight',
                                       'pval_intercept','pval_treatment','pval_status','pval_age','pval_preg','pval_weight')
          
        }
      } else if (LM == "basicPreg+Weight") {
        
        if (singleTrt) {
          
          fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                       fixed_cov$age +
                                       fixed_cov$preg +
                                       fixed_cov$weight)
          
          emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
          rownames(emmremlOutput) <- rownames(resids)
          colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg','beta_weight',
                                       'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg','var_beta_weight',
                                       'pval_intercept','pval_status','pval_age','pval_preg','pval_weight')
        } else {
          fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                       fixed_cov$status + 
                                       fixed_cov$age +
                                       fixed_cov$preg +
                                       fixed_cov$weight)
          
          emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
          rownames(emmremlOutput) <- rownames(resids)
          colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_status','beta_age','beta_preg','beta_weight',
                                       'var_beta_intercept','var_beta_treatment','var_beta_status','var_beta_age','var_beta_preg','var_beta_weight',
                                       'pval_intercept','pval_treatment','pval_status','pval_age','pval_preg','pval_weight')
          
        }
      } else if (LM == "nestedWeight") {
        
        #center the weights by status first
        ###WORKTODO###
        
        
        if (unique(fixed_cov$sex) == "F" ) {
          if (singleTrt) {
            
            fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                         fixed_cov$age +
                                         fixed_cov$preg +
                                         fixed_cov$weightStatusScaled:fixed_cov$status)
            
            emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
            rownames(emmremlOutput) <- rownames(resids)
            colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg','beta_weight:sub','beta_weight:dom',
                                         'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg','var_beta_weight:sub','var_beta_weight:dom',
                                         'pval_intercept','pval_status','pval_age','pval_preg','pval_weight:sub','pval_weight:dom')
    

          } else {
            fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                         fixed_cov$status + 
                                         fixed_cov$age +
                                         fixed_cov$preg +
                                         fixed_cov$weightStatusScaled:fixed_cov$status)
            
            emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
            rownames(emmremlOutput) <- rownames(resids)
            colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_status','beta_age','beta_preg','beta_weight:sub','beta_weight:dom',
                                         'var_beta_intercept','var_beta_treatment','var_beta_status','var_beta_age','var_beta_preg','var_beta_weight:sub','var_beta_weight:dom',
                                         'pval_intercept','pval_treatment','pval_status','pval_age','pval_preg','pval_weight:sub','pval_weight:dom')
            
          }
        } else {
          #males only, drop preg
          if (singleTrt) {
            
            fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                         fixed_cov$age +
                                         fixed_cov$weightStatusScaled:fixed_cov$status)
            
            emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
            rownames(emmremlOutput) <- rownames(resids)
            colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_weight:sub','beta_weight:dom',
                                         'var_beta_intercept','var_beta_status','var_beta_age','var_beta_weight:sub','var_beta_weight:dom',
                                         'pval_intercept','pval_status','pval_age','pval_weight:sub','pval_weight:dom')
          } else {
            fixedCovMM <- model.matrix(~ fixed_cov$treatment + 
                                         fixed_cov$status + 
                                         fixed_cov$age +
                                         fixed_cov$weightStatusScaled:fixed_cov$status)
            
            emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
            rownames(emmremlOutput) <- rownames(resids)
            colnames(emmremlOutput) <- c('beta_intercept','beta_treatment','beta_status','beta_age','beta_weight:sub','beta_weight:dom',
                                         'var_beta_intercept','var_beta_treatment','var_beta_status','var_beta_age','var_beta_weight:sub','var_beta_weight:dom',
                                         'pval_intercept','pval_treatment','pval_status','pval_age','pval_weight:sub','pval_weight:dom')
            
          }
        }
      } else if (LM == "basicSeason") {
        
        fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                     fixed_cov$age +
                                     fixed_cov$season)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_season',
                                     'var_beta_intercept','var_beta_status','var_beta_age','var_beta_season',
                                     'pval_intercept','pval_status','pval_age','pval_season')
      } else if (LM == "allTrtF") {
        fixedCovMM <- model.matrix(~ fixed_cov$status +
                                     fixed_cov$age +
                                     fixed_cov$treatment)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_trtDex','beta_trtGard','beta_trtLPS',
                                     'var_beta_intercept','var_beta_status','var_beta_age','var_beta_trtDex','var_beta_trtGard','var_beta_trtLPS',
                                     'pval_intercept','pval_status','pval_age','pval_trtDex','pval_trtGard','pval_trtLPS')
      } else if (LM == "allTrtM") {
        fixedCovMM <- model.matrix(~ fixed_cov$status +
                                     fixed_cov$age +
                                     fixed_cov$treatment)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_trtDex','beta_trtGard','beta_trtLPS',
                                     'var_beta_intercept','var_beta_status','var_beta_age','var_beta_trtDex','var_beta_trtGard','var_beta_trtLPS',
                                     'pval_intercept','pval_status','pval_age','pval_trtDex','pval_trtGard','pval_trtLPS')
      } else if (LM == "allTrtF+Weight") {
        fixedCovMM <- model.matrix(~ fixed_cov$status +
                                     fixed_cov$age +
                                     fixed_cov$weight +
                                     fixed_cov$treatment)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_weight','beta_trtDex','beta_trtGard','beta_trtLPS',
                                     'var_beta_intercept','var_beta_status','var_beta_age','var_beta_weight','var_beta_trtDex','var_beta_trtGard','var_beta_trtLPS',
                                     'pval_intercept','pval_status','pval_age','pval_weight','pval_trtDex','pval_trtGard','pval_trtLPS')
      } else if (LM == "allTrtF+Preg") {
        fixedCovMM <- model.matrix(~ fixed_cov$status +
                                     fixed_cov$age +
                                     fixed_cov$preg +
                                     fixed_cov$treatment)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg','beta_trtDex','beta_trtGard','beta_trtLPS',
                                     'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg','var_beta_trtDex','var_beta_trtGard','var_beta_trtLPS',
                                     'pval_intercept','pval_status','pval_age','pval_preg','pval_trtDex','pval_trtGard','pval_trtLPS')
      } else if (LM == "allTrtF+Preg+Weight") {
        fixedCovMM <- model.matrix(~ fixed_cov$status +
                                     fixed_cov$age +
                                     fixed_cov$preg +
                                     fixed_cov$treatment +
                                     fixed_cov$weight)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg','beta_trtDex','beta_trtGard','beta_trtLPS','beta_weight',
                                     'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg','var_beta_trtDex','var_beta_trtGard','var_beta_trtLPS','var_beta_weight',
                                     'pval_intercept','pval_status','pval_age','pval_preg','pval_trtDex','pval_trtGard','pval_trtLPS','pval_weight')
      } else if (LM == "allTrtM+Weight") {
        fixedCovMM <- model.matrix(~ fixed_cov$status +
                                     fixed_cov$age +
                                     fixed_cov$treatment +
                                     fixed_cov$weight)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_trtDex','beta_trtGard','beta_trtLPS','beta_weight',
                                     'var_beta_intercept','var_beta_status','var_beta_age','var_beta_trtDex','var_beta_trtGard','var_beta_trtLPS','var_beta_weight',
                                     'pval_intercept','pval_status','pval_age','pval_trtDex','pval_trtGard','pval_trtLPS','pval_weight')
      } else if (LM == "allTrtMF+Weight") {
        fixedCovMM <- model.matrix(~ fixed_cov$status +
                                     fixed_cov$age +
                                     fixed_cov$treatment +
                                     fixed_cov$weight +
                                     fixed_cov$sex)
        
        emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
        rownames(emmremlOutput) <- rownames(resids)
        colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_trtDex','beta_trtGard','beta_trtLPS','beta_weight','beta_sex',
                                     'var_beta_intercept','var_beta_status','var_beta_age','var_beta_trtDex','var_beta_trtGard','var_beta_trtLPS','var_beta_weight','var_beta_sex',
                                     'pval_intercept','pval_status','pval_age','pval_trtDex','pval_trtGard','pval_trtLPS','pval_weight','pval_sex')
      }
  
  #vUcheck <- vector()
  #vEcheck <- vector()
  
  for (g in 1:nGenes) {
      temp <- emmreml(y = resids[g,],
                      X = fixedCovMM,
                      K = identityK,
                      Z = sampleIndMap,
                      varbetahat = T,varuhat=T,PEVuhat=T,test=T)
      
      p <- temp$pvalbeta[,"none"]
      varb <- temp$varbetahat
      b <- temp$betahat
      
      #vUcheck <- c(vUcheck,temp$Vu)
      #vEcheck <- c(vEcheck,temp$Ve)
      
      emmremlOutput[g,] <- as.vector(c(b,varb,p))
  }
  #system("say dingh")
  
  #ggplot() + geom_histogram(aes(x = emmremlOutput$pval_status), bins = 100)
  #ggplot() + geom_histogram(aes(x = emmremlOutput$pval_weight), bins = 100)
  #ggplot() + geom_point(aes(x = vUcheck, y = vEcheck))
  #sum(vUcheck > .5)
  #ggplot() + geom_point(aes(x = emmremlOutput$beta_weight, y = emmremlOutput$beta_status), alpha = .2)
  #cor(emmremlOutput$beta_status, emmremlOutput$beta_weight)
  #vUcheckPerm2rand <- vUcheck
  #vEcheckPerm2rand <- vEcheck
  
  #ggplot() + geom_point(aes(x = vUcheckPerm0, y = vUcheckPerm2rand), alpha = .2) + geom_abline(slope = 1, intercept = 0)
  #ggplot() + geom_point(aes(x = vEcheckPerm0, y = vEcheckPerm2rand), alpha = .2) + geom_abline(slope = 1, intercept = 0)
  
  #add some model info to the titles
  if (trtName %in% c("LPS","Gard","Dex")) {
    if (singleTrt) {
      trtTitle <- paste0(trtName,"only")
    } else {
      trtTitle <- paste0(trtName,"-Ctrl")
    }
  } else {
    trtTitle <- "Ctrlonly"
  }
  
  #if using a Kmat note it in the model
  if (Kmat) {
    LMbasic <- LM
    LM <- paste0(LM,"wKmatrix")
  } else {
    LMbasic <- LM
  }
  
    
  #figOut <- paste0(figDir,modelName)
  if ( perms == FALSE ) {
    ggplot() +
      geom_histogram(aes(x = emmremlOutput$pval_intercept), bins = 100) +
      xlab("p-value - intercept") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
    ggsave(paste0(modelDir,"/",LM,"_pval_intercept.png"), width = 6, height = 4)
    
    if (singleTrt == FALSE & combineTrt == FALSE) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_treatment), bins = 100) +
        xlab("p-value - treatment") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_treatment.png"), width = 6, height = 4)
    } else if (combineTrt == TRUE) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_trtDex), bins = 100) +
        xlab("p-value - Dex trt") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_trtDex.png"), width = 6, height = 4)
      
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_trtGard), bins = 100) +
        xlab("p-value - Gard trt") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_trtGard.png"), width = 6, height = 4)
      
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_trtLPS), bins = 100) +
        xlab("p-value - LPS trt") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_trtLPS.png"), width = 6, height = 4)
    }
    
    ggplot() +
      geom_histogram(aes(x = emmremlOutput$pval_age), bins = 100) +
      xlab("p-value - age") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
    ggsave(paste0(modelDir,"/",LM,"_pval_age.png"), width = 6, height = 4)
    
    if (LMbasic %in% c("basic","basic+Preg","nestedPreg","allTrtF","allTrtF+Weight","allTrtM",
                       "allTrtM+Weight","basic+Preg+Weight","basic+Weight","nestedWeight","allTrtF+Preg+Weight")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_status), bins = 100) +
        xlab("p-value - status") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_status.png"), width = 6, height = 4)
    }
    
    if (LMbasic %in% c("basic+Preg","nestedStatus","nestedWPreg","nested+Preg+Weight","allTrtF+Preg",
                       "nestedWeight","basic+Preg+Weight","basicPregWeightOnly","allTrtF+Preg+Weight")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_preg), bins = 100) +
        xlab("p-value - preg") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_preg.png"), width = 6, height = 4)
      
      #ggplot() +
      #  geom_histogram(aes(x = emmremlOutput$pval_lact), bins = 100) +
      #  xlab("p-value - lact") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      #ggsave(paste0(modelDir,"/",LM,"_pval_lact.png"), width = 6, height = 4)
    }
    
    if (LMbasic %in% c("nestedPreg")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_Dom.Preg), bins = 100) +
        xlab("p-value - Dom:Preg") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_DomxPreg.png"), width = 6, height = 4)
      
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_Sub.Preg), bins = 100) +
        xlab("p-value - Sub:Preg") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_SubxPreg.png"), width = 6, height = 4)
    }
    
    if (LMbasic %in% c("nestedStatus")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_Preg.Status), bins = 100) +
        xlab("p-value - Preg:Status") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_DomxPreg.png"), width = 6, height = 4)
      
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_NonPreg.Status), bins = 100) +
        xlab("p-value - NotPreg:Status") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_SubxPreg.png"), width = 6, height = 4)
    }
    
    if (LMbasic %in% c("basicSeason")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_season), bins = 100) +
        xlab("p-value - season") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_season.png"), width = 6, height = 4)
      
    }
      
    if (LMbasic %in% c("nested","nested+Weight","nestedWPreg","nested+Preg+Weight")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$`pval_Ctrl:Status`), bins = 100) +
        xlab("p-value - Ctrl:Status") + 
        ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_CtrlxStatus.png"), width = 6, height = 4)
      
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$`pval_Trt:Status`), bins = 100) +
        xlab("p-value - Trt:Status") + 
        ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_TrtxStatus.png"), width = 6, height = 4)
    }
    
    #weight pvals
    if (LMbasic %in% c("nested+Weight","basic+Preg+Weight","nested+Preg+Weight","basic+Weight",
                       "allTrtF+Weight","allTrtM+Weight","basicPregWeightOnly","allTrtF+Preg+Weight")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_weight), bins = 100) +
        xlab("p-value - weight") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_weight.png"), width = 6, height = 4)
      
    }
    
    if (LMbasic %in% c("nestedWeight")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$`pval_weight:dom`), bins = 100) +
        xlab("p-value - Dom:Weight") + 
        ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_DomxWeight.png"), width = 6, height = 4)
      
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$`pval_weight:sub`), bins = 100) +
        xlab("p-value - Sub:Weight") + 
        ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_SubxWeight.png"), width = 6, height = 4)
    }
    
    
    #write out the emmreml
    write.table(x = emmremlOutput, file = paste0(outDir,"/",gatherSubset,"_model",model,"_",LM,"_emmreml"))
    
  } else {
    if ( permVar == "treatment" ) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_treatment), bins = 100) +
        xlab("p-value - treatment") + 
        ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName," | perm",formatC(permNum, width=3, flag="0")))
      
      ggsave(paste0(modelDir,"/",LM,"_pval_treatment_perm",
                    formatC(permNum, width=3, flag="0"),".png"), width = 6, height = 4)
    } else if ( permVar == "status" ) {
      if (LMbasic == "basic") {
        ggplot() +
          geom_histogram(aes(x = emmremlOutput$pval_status), bins = 100) +
          xlab("p-value - status") + 
          ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName," | perm",formatC(permNum, width=3, flag="0")))
        
        ggsave(paste0(modelDir,"/",LM,"_pval_status_perm",
                      formatC(permNum, width=3, flag="0"),".png"), width = 6, height = 4)
      }
      if (LMbasic %in% c("nested","nested+Weight")) {
        ggplot() +
          geom_histogram(aes(x = emmremlOutput$pval_Trt.Status), bins = 100) +
          xlab("p-value - Trt:Status") + 
          ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName," | perm",formatC(permNum, width=3, flag="0")))
        
        ggsave(paste0(modelDir,"/",LM,"_pval_TrtxStatus_perm",
                      formatC(permNum, width=3, flag="0"),".png"), width = 6, height = 4)
        ggplot() +
          geom_histogram(aes(x = emmremlOutput$pval_Ctrl.Status), bins = 100) +
          xlab("p-value - Ctrl:Status") + 
          ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName," | perm",formatC(permNum, width=3, flag="0")))
        
        ggsave(paste0(modelDir,"/",LM,"_pval_CtrlxStatus_perm",
                      formatC(permNum, width=3, flag="0"),".png"), width = 6, height = 4)
        
      }
      
    }
    
    write.table(x = emmremlOutput, file = paste0(outDir,"/",gatherSubset,"_model",model,"_",LM,"_emmreml_",
                                                 permVar,"_v",permVersion,"_perm",formatC(permNum, width=3, flag="0")))
    ### saved permuted fixed_cov
    saveRDS(object = fixed_cov, file = paste0(outDir,"/",gatherSubset,"_metadata_",permVar,"_v",permVersion,"_perm",
                                              formatC(permNum, width=3, flag="0"),".RDS"))
  }
}
