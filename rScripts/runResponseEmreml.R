## using counts and metadata, generate residuals for emmreml analysis

# test a change
# 
batchList <- c("allBatches")
focalSex <- "F"
focalTreatment <- "LPS"
minGenes <- 10000
workDir <- "~/Dropbox (Personal)/meerkats/RNA/meerkatGE/"
minLogCPM <- 5
# 
Kmat <- TRUE
Kfile <- "~/Dropbox (Personal)/meerkats/KMPsamples/KMPrelatedness_matrix.RDS"
# 
# # 
gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"_",minGenes,"genes_",minLogCPM,"lCPM")
# # 
model <- "F"
LM <- "basicResp+Preg+Weight"
# LM <- "basicResp"
# perms <- FALSE

# perms <- TRUE
# permVar <- "status"
# permNum <- 7



runResponseEmreml <- function(workDir, gatherSubset, model, LM, perms = FALSE, permVar = "treatment", permNum = 1,
                              Kmat = FALSE, Kfile = "~/Dropbox (Personal)/meerkats/KMPsamples/KMPrelatedness_matrix.RDS"){
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
  resids <- readRDS(file = paste0(outDir,"/",gatherSubset,"_model",model,"_resids.RDS"))
  
  if (Kmat) {
    K <- readRDS(Kfile)
  }
  
  #determine focal treatment
  focalTreatment <- subset(subsetMeta, !Trt == "Ctrl")$Trt[1]
  
  responseCaptures <- names(table(subsetMeta$captureref)[table(subsetMeta$captureref) > 1])
  
  responseCaptBoth <- vector()
  
  #r <- responseCaptures[14]
  for ( r in responseCaptures ) {
    capTrts <- subset(subsetMeta, captureref == r)$Trt
    if ( "Ctrl" %in% capTrts & focalTreatment %in% capTrts ) {
      responseCaptBoth <- c(responseCaptBoth,r)
    }
  }
  
  #use response captures to narrow down meta & residuals
  respMeta <- subset(subsetMeta, captureref %in% responseCaptBoth)
  respResids <- resids[,colnames(resids) %in% respMeta$sampleID]
  
  #split both into LPS & ctrl files
  all(colnames(respResids) == respMeta$sampleID)
  
  #drop samples with repeats, for simplicity
  if (!all(table(respMeta$captureref) == 2)) {
    dropCaptures <- names(table(respMeta$captureref)[table(respMeta$captureref) > 2])
    respMeta <- subset(respMeta, !captureref %in% dropCaptures)
    respResids <- resids[,colnames(resids) %in% respMeta$sampleID]
  }
  
  #subtract LPS - Ctrl
  trtMeta <- subset(respMeta, Trt == focalTreatment)
  ctrlMeta <- subset(respMeta, Trt == "Ctrl")
  
  all(trtMeta$captureref == ctrlMeta$captureref)
  
  ctrlResids <- respResids[,colnames(respResids) %in% ctrlMeta$sampleID]
  trtResids <- respResids[,colnames(respResids) %in% trtMeta$sampleID]
  
  respResids <- trtResids - ctrlResids
  ctrlMeta$Trt <- "Response"
  ctrlMeta$treatment.x <- "resp"
  
  respMeta <- ctrlMeta
  subsetMeta <- respMeta
  
  #save data
  saveRDS(object = respMeta, file = paste0(outDir,"/",gatherSubset,"_responseMetadata.RDS"))
  saveRDS(object = respResids, file = paste0(outDir,"/",gatherSubset,"_model",model,"_responseResids.RDS"))
  resids <- respResids
  
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
  
  fixed_cov <- subsetMeta[which(colnames(subsetMeta) %in% c("sex","status_dom_sub","ageAtSample","permID","animal_id","pregLact","climSeason","avgWeight6mo"))]
  
  colnames(fixed_cov)[which(colnames(fixed_cov) == "Trt")] <- "treatment"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "status_dom_sub")] <- "status"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "ageAtSample")] <- "age"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "pregLact")] <- "reprod"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "climSeason")] <- "season"
  colnames(fixed_cov)[which(colnames(fixed_cov) == "avgWeight6mo")] <- "weight"
  
  trtName <- focalTreatment

  fixed_cov$age <- scale(fixed_cov$age, scale = F, center = T)
  fixed_cov$weightUnScaled <- fixed_cov$weight
  animalWeights <- unique(fixed_cov$weight)
  scaledWeights <- scale(animalWeights, scale = F, center = T)
  
  for (n in 1:dim(fixed_cov)[1]) {
    matchRow <- which(animalWeights == fixed_cov$weightUnScaled[n])
    fixed_cov$weight[n] <- scaledWeights[matchRow]
  }
  
  fixed_cov$status <- factor(fixed_cov$status, levels = c("S","D"))
  
  fixed_cov$reprod <- factor(fixed_cov$reprod, levels = c("none","preg","lact"))
  fixed_cov$preg <- ifelse(fixed_cov$reprod == "preg", TRUE, FALSE)
  
  fixed_cov$season <- factor(fixed_cov$season, levels = c("warm-wet","cold-dry"))
  
  #### need some extra help here to identify individuals in the fixed_cov object
  if ( perms ) {
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
      
    } else if ( permVar == "status" ) {
      fixed_cov$animalStatus <- paste0(fixed_cov$animal_id,"_",fixed_cov$status)
      
      #permute status variable here      
      possibleStatus <- unique(fixed_cov$status)
      
      anmlStat <- unique(fixed_cov$animalStatus)[14]
      for (anmlStat in unique(fixed_cov$animalStatus)) {
        covRows <- which(fixed_cov$animalStatus == anmlStat)
        
        actualStatus <- fixed_cov$status[covRows]
        #permTrts <- actualTrts
        
        #randomly flip or not
        if (runif(1) < .5) {
          permStatus <- rep(possibleStatus[1],length(actualStatus))
        } else {
          permStatus <- rep(possibleStatus[2],length(actualStatus))
        }
        
        fixed_cov$status[covRows] <- permStatus
      }
    } else if ( permVar == "preg") {
      ### PERMUTE PREGNANCY STATUS FOR EACH ANIMAL
      #permute treatment variable here
      
      ## this is pasted in from treatment... but needs to be looked at
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
    }
  }
  
  if ( LM == "basicResp" ) {
    
      #basic model
      fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                   fixed_cov$age)
      
      emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
      rownames(emmremlOutput) <- rownames(resids)
      colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age',
                                   'var_beta_intercept','var_beta_status','var_beta_age',
                                   'pval_intercept','pval_status','pval_age')

  } else if ( LM == "basicResp+Weight" ) {
    
      #basic model
      fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                   fixed_cov$age +
                                   fixed_cov$weight)
      
      emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
      rownames(emmremlOutput) <- rownames(resids)
      colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_weight',
                                   'var_beta_intercept','var_beta_status','var_beta_age','var_beta_weight',
                                   'pval_intercept','pval_status','pval_age','pval_weight')
    
  } else if (LM == "basicResp+Preg") {
    
    fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                 fixed_cov$age + 
                                 fixed_cov$preg)
    
    emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
    rownames(emmremlOutput) <- rownames(resids)
    colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg',
                                 'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg',
                                 'pval_intercept','pval_status','pval_age','pval_preg')
    
  } else if (LM == "basicResp+Preg+Weight") {
    
    fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                 fixed_cov$age +
                                 fixed_cov$preg +
                                 fixed_cov$weight)
    
    emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
    rownames(emmremlOutput) <- rownames(resids)
    colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg','beta_weight',
                                 'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg','var_beta_weight',
                                 'pval_intercept','pval_status','pval_age','pval_preg','pval_weight')
  } else if (LM == "respNestedPreg") {
    
    fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                 fixed_cov$age  +
                                 fixed_cov$weight +
                                 fixed_cov$preg:fixed_cov$status)
    
    emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
    rownames(emmremlOutput) <- rownames(resids)
    colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_weight','beta_Sub:Preg','beta_Dom:Preg',
                                 'var_beta_intercept','var_beta_status','var_beta_age','var_beta_weight','var_beta_Sub:Preg','var_beta_Dom:Preg',
                                 'pval_intercept','pval_status','pval_age','pval_weight','pval_Sub:Preg','pval_Dom:Preg')
  } else if (LM == "respNestedWeight") {
    
    fixedCovMM <- model.matrix(~ fixed_cov$status + 
                                 fixed_cov$age  +
                                 fixed_cov$preg +
                                 fixed_cov$weight:fixed_cov$status)
    
    emmremlOutput <- as.data.frame(diag(0,ncol = dim(fixedCovMM)[2] * 3, nrow = dim(resids)[1]))
    rownames(emmremlOutput) <- rownames(resids)
    colnames(emmremlOutput) <- c('beta_intercept','beta_status','beta_age','beta_preg','beta_Sub:Weight','beta_Dom:Weight',
                                 'var_beta_intercept','var_beta_status','var_beta_age','var_beta_preg','var_beta_Sub:Weight','var_beta_Dom:Weight',
                                 'pval_intercept','pval_status','pval_age','pval_preg','pval_Sub:Weight','pval_Dom:Weight')
  }
  
  for (g in 1:nGenes) {
    temp <- emmreml(y = resids[g,],
                    X = fixedCovMM,
                    K = identityK,
                    Z = sampleIndMap,
                    varbetahat = T,varuhat=T,PEVuhat=T,test=T)
    
    p <- temp$pvalbeta[,"none"]
    varb <- temp$varbetahat
    b <- temp$betahat
    
    emmremlOutput[g,] <- as.vector(c(b,varb,p))
  }
  
  #add some model info to the titles
  trtTitle <- paste0(trtName,"-Response")
  
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
    
    ggplot() +
      geom_histogram(aes(x = emmremlOutput$pval_age), bins = 100) +
      xlab("p-value - age") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
    ggsave(paste0(modelDir,"/",LM,"_pval_age.png"), width = 6, height = 4)
    
    ggplot() +
      geom_histogram(aes(x = emmremlOutput$pval_status), bins = 100) +
      xlab("p-value - status") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
    ggsave(paste0(modelDir,"/",LM,"_pval_status.png"), width = 6, height = 4)
    
    if (LM %in% c("basicResp+Preg","basicResp+Preg+Weight","respNestedWeight")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_preg), bins = 100) +
        xlab("p-value - preg") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_preg.png"), width = 6, height = 4)
    }
    
    if (LM %in% c("basicResp+Weight","basicResp+Preg+Weight","nestedResp+PregStatus")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$pval_weight), bins = 100) +
        xlab("p-value - weight") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_weight.png"), width = 6, height = 4)
    }
    
    if (LM %in% c("respNestedPreg")) {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$`pval_Dom:Preg`), bins = 100) +
        xlab("p-value - Dom:Preg") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_DomxPreg.png"), width = 6, height = 4)
      
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$`pval_Sub:Preg`), bins = 100) +
        xlab("p-value - Sub:Preg") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_SubxPreg.png"), width = 6, height = 4)
    }
    
    if (LM == "respNestedWeight") {
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$`pval_Dom:Weight`), bins = 100) +
        xlab("p-value - Dom:Weight") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_DomxWeight.png"), width = 6, height = 4)
      
      ggplot() +
        geom_histogram(aes(x = emmremlOutput$`pval_Sub:Weight`), bins = 100) +
        xlab("p-value - Sub:Weight") + ggtitle(paste0("resid - ",model," | LM - ",LM," | ",trtName))
      ggsave(paste0(modelDir,"/",LM,"_pval_SubxWeight.png"), width = 6, height = 4)
    }
    
    write.table(x = emmremlOutput, file = paste0(outDir,"/",gatherSubset,"_model",model,"_",LM,"_emmreml"))
    
    } else {
      #perm plots
      
      
      #perm output
      write.table(x = emmremlOutput, file = paste0(outDir,"/",gatherSubset,"_model",model,"_",LM,"_emmreml_",
                                                   permVar,"_perm",formatC(permNum, width=3, flag="0")))
      ### saved permuted fixed_cov
      saveRDS(object = fixed_cov, file = paste0(outDir,"/",gatherSubset,"_metadata_",permVar,"_perm",
                                                formatC(permNum, width=3, flag="0"),".RDS"))
      
    }
}
