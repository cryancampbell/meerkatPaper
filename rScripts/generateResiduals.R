## 
## ---------------------------
##
## Script name: generateResiduals.R
##
## Purpose of script: using counts and metadata, generate residuals for emmreml analysis
##
## Author: CRC
##
## Date Created: 2022-02-03
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

## set working directory for Mac and PC

#?#?#
# how do libraries work? shoudl i load them in here as well?

# workDir <- workDir
# batchList <- c("allBatches")
# focalSex <- "M"
# focalTreatment <- "LPS"
# minLogCPM <- 5
# minGenes <- 10000
# singleTrt <- FALSE
# combineTrt <- FALSE
#  
#  if (singleTrt) {
#    gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"only_",minGenes,"genes_",minLogCPM,"lCPM")
#  } else if (combineTrt) {
#    gatherSubset <- paste0(batchList,"_",focalSex,"_allTrt_",minGenes,"genes_",minLogCPM,"lCPM")
#  } else {
#    gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"_",minGenes,"genes_",minLogCPM,"lCPM")
#  }
#  
#  model <- "F"

#A: model i used with 2020 data, correct for Total Reads, Unique Genes, Qubit Conc, Group
#B: 2021 data, correct for Total Reads, Unique Genes, Percent of Reads In Feature, Percent Mapped to Genome
#C: B + Library Batch data
#D: B + Library Batch + Qubit Conc
#E: B + Seq Batch (combined only)
#F: B + Library Batch - Reads
#G: B + Library Batch + Qubit - Reads OR D - Reads

#K: Combat seq batches
#L: Combat seq batch + technical correction
#M: Combat library batches
#N: Combat library batches + technical correction

#P: F + Preg
#Q: G + Preg

#R: F + Status
#S: G + Status

#use G, Q, & S? for "allTrt"

#T: F + Weight
#U: G + Weight

#V: F + Status Centered Weight
#W: G + Status Centered Weight

#H: F + Age
#I: G + Age

#X: F + Status-Centered Age
#Y: G + Status-Centered Age


workDir <- workDir
batchList <- c("allBatches")
focalSex <- "F"
focalTreatment <- "LPS"
minLogCPM <- 5
minGenes <- 10000
#if TRUE, focal treatment is the ONLY treatment used
#focal treatment object becomes "LPSonly" or "Ctrlonly" instead of "LPS"
singleTrt <- TRUE
#combineTrt is all 4 treatments (LPS/Gard/Dex/Ctrl) analyzed jointly
combineTrt <- FALSE

### longitudinal
longReSample <- TRUE
ReSampVar <- "status"
model <- "F"

gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"only_",minGenes,"genes_",minLogCPM,"lCPM_",ReSampVar,"ReSamp")


genResids <- function(workDir, gatherSubset, model){

  dataDir <- paste0(workDir,"data/")
  outDir <- paste0(workDir,"output/")
  figDir <- paste0(workDir,"figures/")
  
  #read in data
  subsetMeta <- readRDS(file = paste0(outDir,"/",gatherSubset,"_metadata.RDS"))
  subsetCountsDF <- readRDS(file = paste0(outDir,"/",gatherSubset,"_counts.RDS"))
  filterCountsDF <- readRDS(file = paste0(outDir,"/",gatherSubset,"_filteredCounts.RDS"))
  
  for (lmMethod in c("none","used")) {
  #run PCA without any correction first, always:
    if (lmMethod == "none") {
      lmUsed <- "none"
      titleText <- " (directly on GE)"
      
      #normalize the counts
      vNorm <- voom(filterCountsDF,plot=FALSE)
      pca <- prcomp(cor(vNorm$E, method = "spearman"))
      
      df_out <- as.data.frame(pca$x)
      dfMeta <- cbind(df_out,subsetMeta)
      pcaNoLM <- pca
      summary(pca)
    } else {
      #then run a correction of some type, model A-...

      ## if the model DOES include combat, run that first: guesses at letters
      
      if ( model %in% c("K","L","M","N") ) {
        if (model %in% c("K","L")) {
          #RUN COMBAT
          lmUsed <- "ComBat"
          
          titleText <- " (ComBat Used)"
          
          vNorm <- voom(filterCountsDF,plot=FALSE)
          
          expressionMatrix <- as.matrix(vNorm$E)
          
          ##### THIS IS WHERE WE EMPLOY COMBAT - vNormLM
          #combat
          batch <- subsetMeta$seqBatch
          
          #groupMulti <- subsetMeta[,colnames(subsetMeta) %in% c("Trt","status_dom_sub")]
          groupMulti <- model.matrix(~1 + subsetMeta$Trt + subsetMeta$status_dom_sub)
          
          combat_edata <- ComBat(dat=expressionMatrix,
                                 batch=batch, 
                                 mod=groupMulti, 
                                 par.prior=TRUE, 
                                 prior.plots=FALSE)
          
          if (model == "K") {
            residuals <- combat_edata
          } else if (model == "L") {
            #control for tech variables
            design <- model.matrix(~ subsetMeta$Reads + 
                                     subsetMeta$UniqueGenes +
                                     subsetMeta$PercentInFeat +
                                     subsetMeta$percentMapped +
                                     subsetMeta$qubitConc)
            lmUsed <- "ComBat_wLM"
            
            titleText <- " (ComBat + lm)"
            
            #no weights
            
            #### then fit the design to the combat output
            fit <- lmFit(combat_edata,design)
            residuals <- residuals.MArrayLM(object=fit, combat_edata)
          }
          
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          pcaNoLM <- pca
          summary(pca)
          
        } else if (model %in% c("M","N")) {
              #RUN COMBAT
              lmUsed <- "ComBat"
              
              titleText <- " (ComBat Used)"
              
              vNorm <- voom(filterCountsDF,plot=FALSE)
              
              expressionMatrix <- as.matrix(vNorm$E)
          
              ##### THIS IS WHERE WE EMPLOY COMBAT - vNormLM
              #combat
              batch <- subsetMeta$seqLibBatch
              
              #groupMulti <- subsetMeta[,colnames(subsetMeta) %in% c("Trt","status_dom_sub")]
              groupMulti <- model.matrix(~1 + subsetMeta$Trt + subsetMeta$status_dom_sub)
              
              combat_edata <- ComBat(dat=expressionMatrix,
                                     batch=batch, 
                                     mod=groupMulti, 
                                     par.prior=TRUE, 
                                     prior.plots=FALSE)
              
              if (model == "M") {
                residuals <- combat_edata
              } else if (model == "N") {
                #control for tech variables
                design <- model.matrix(~ subsetMeta$Reads + 
                                         subsetMeta$UniqueGenes +
                                         subsetMeta$PercentInFeat +
                                         subsetMeta$percentMapped +
                                         subsetMeta$qubitConc)
                lmUsed <- "ComBat_wLM"
                
                titleText <- " (ComBat + lm)"
                
                #no weights
                
                #### then fit the design to the combat output
                fit <- lmFit(combat_edata,design)
                residuals <- residuals.MArrayLM(object=fit, combat_edata)
              }
              
              
              pca <- prcomp(cor(residuals, method = "spearman"))
              
              df_out <- as.data.frame(pca$x)
              dfMeta <- cbind(df_out,subsetMeta)
              pcaNoLM <- pca
              summary(pca)
            
            }
            
          }
          
      
  
  ## if the model does NOT include combat, start here with a typical PCA with / without LM
      if ( model %in% c("A","B","C","D","E","F","G","P","Q","R","S","T","U","V","W","H","I","X","Y") ) {
        #control for tech variables
        # find it by model
        if (model == "A") {
          design <- model.matrix(~ subsetMeta$Reads + subsetMeta$UniqueGenes +
                                   subsetMeta$qubitConc + subsetMeta$group)
          
          #from 2020 model
          #metaSubset$inputReads + metaSubset$genesMapped + 
          #metaSubset$qubitConc + metaSubset$group
          
          lmUsed <- "Reads_Genes_Qubit_Group"

          titleText <- " (lm Used)"
        
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
        
          pcaWithLM <- pca
          summary(pca)
          
        } else if (model == "B") {
          design <- model.matrix(~ subsetMeta$Reads + subsetMeta$UniqueGenes +
                                   subsetMeta$PercentInFeat + subsetMeta$percentMapped)
          
          #from 2021 model
          #modelMeta$uniqueGenes + modelMeta$Reads +
          #modelMeta$percentInFeat + modelMeta$percentMapped

          lmUsed <- "Reads_Genes_PercentInFeat_PercentMapped"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "C") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$Reads +
                                   subsetMeta$PercentInFeat + subsetMeta$percentMapped + subsetMeta$seqLibBatch)
          
          lmUsed <- "Reads_Genes_PercentInFeat_PercentMapped_libBatch"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "D") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$Reads + subsetMeta$qubitConc +
                                   subsetMeta$PercentInFeat + subsetMeta$percentMapped + subsetMeta$seqLibBatch)
          
          lmUsed <- "Reads_Genes_PercentInFeat_PercentMapped_qubit_libBatch"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "E") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$Reads +
                                   subsetMeta$PercentInFeat + subsetMeta$percentMapped + subsetMeta$seqBatch)
          
          lmUsed <- "Reads_Genes_PercentInFeat_PercentMapped_seqBatch"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "F") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + subsetMeta$seqLibBatch)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "G") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                   subsetMeta$seqLibBatch + subsetMeta$qubitConc)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_qubit"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "P") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                   subsetMeta$seqLibBatch + subsetMeta$pregStatusBackDate35)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_Preg"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "Q") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                   subsetMeta$seqLibBatch + subsetMeta$qubitConc + subsetMeta$pregStatusBackDate35)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_qubit_Preg"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "R") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                   subsetMeta$seqLibBatch + subsetMeta$status_dom_sub)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_Status"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "S") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                   subsetMeta$seqLibBatch + subsetMeta$status_dom_sub)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_Status"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "T") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                   subsetMeta$seqLibBatch + subsetMeta$centWeight)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_qubit_Weight"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
        } else if (model == "U") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                   subsetMeta$seqLibBatch + subsetMeta$qubitConc + subsetMeta$centWeight)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_Weight"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
          } else if (model == "V") {
          design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                   subsetMeta$seqLibBatch + subsetMeta$statCentWeight)
          
          lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_qubit_centWeight"
          
          titleText <- " (lm Used)"
          
          vNormLM <- voom(filterCountsDF,design,plot=FALSE)
          #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
          
          fit <- lmFit(vNormLM,design)
          residuals <- residuals.MArrayLM(object=fit, vNormLM)
          
          pca <- prcomp(cor(residuals, method = "spearman"))
          
          df_out <- as.data.frame(pca$x)
          dfMeta <- cbind(df_out,subsetMeta)
          
          pcaWithLM <- pca
          summary(pca)
          } else if (model == "W") {
            design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                     subsetMeta$seqLibBatch + subsetMeta$qubitConc + subsetMeta$statCentWeight)
            
            lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_qubit_centWeight"
            
            titleText <- " (lm Used)"
            
            vNormLM <- voom(filterCountsDF,design,plot=FALSE)
            #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
            
            fit <- lmFit(vNormLM,design)
            residuals <- residuals.MArrayLM(object=fit, vNormLM)
            
            pca <- prcomp(cor(residuals, method = "spearman"))
            
            df_out <- as.data.frame(pca$x)
            dfMeta <- cbind(df_out,subsetMeta)
            
            pcaWithLM <- pca
            summary(pca)
          } else if (model == "H") {
            design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                     subsetMeta$seqLibBatch + subsetMeta$centAge)
            
            lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_age"
            
            titleText <- " (lm Used)"
            
            vNormLM <- voom(filterCountsDF,design,plot=FALSE)
            #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
            
            fit <- lmFit(vNormLM,design)
            residuals <- residuals.MArrayLM(object=fit, vNormLM)
            
            pca <- prcomp(cor(residuals, method = "spearman"))
            
            df_out <- as.data.frame(pca$x)
            dfMeta <- cbind(df_out,subsetMeta)
            
            pcaWithLM <- pca
            summary(pca)
          } else if (model == "I") {
            design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                     subsetMeta$seqLibBatch + subsetMeta$qubitConc + subsetMeta$centAge)
            
            lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_qubit_age"
            
            titleText <- " (lm Used)"
            
            vNormLM <- voom(filterCountsDF,design,plot=FALSE)
            #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
            
            fit <- lmFit(vNormLM,design)
            residuals <- residuals.MArrayLM(object=fit, vNormLM)
            
            pca <- prcomp(cor(residuals, method = "spearman"))
            
            df_out <- as.data.frame(pca$x)
            dfMeta <- cbind(df_out,subsetMeta)
            
            pcaWithLM <- pca
            summary(pca)
          } else if (model == "X") {
            design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                     subsetMeta$seqLibBatch + subsetMeta$statCentAge)
            
            lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_centAge"
            
            titleText <- " (lm Used)"
            
            vNormLM <- voom(filterCountsDF,design,plot=FALSE)
            #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
            
            fit <- lmFit(vNormLM,design)
            residuals <- residuals.MArrayLM(object=fit, vNormLM)
            
            pca <- prcomp(cor(residuals, method = "spearman"))
            
            df_out <- as.data.frame(pca$x)
            dfMeta <- cbind(df_out,subsetMeta)
            
            pcaWithLM <- pca
            summary(pca)
          } else if (model == "Y") {
            design <- model.matrix(~ subsetMeta$UniqueGenes + subsetMeta$PercentInFeat + subsetMeta$percentMapped + 
                                     subsetMeta$seqLibBatch + subsetMeta$qubitConc + subsetMeta$statCentAge)
            
            lmUsed <- "Genes_PercentInFeat_PercentMapped_libBatch_qubit_centAge"
            
            titleText <- " (lm Used)"
            
            vNormLM <- voom(filterCountsDF,design,plot=FALSE)
            #vNormLM <- voomWithQualityWeights(modelCountsLCPMfilt,design,plot=FALSE)
            
            fit <- lmFit(vNormLM,design)
            residuals <- residuals.MArrayLM(object=fit, vNormLM)
            
            pca <- prcomp(cor(residuals, method = "spearman"))
            
            df_out <- as.data.frame(pca$x)
            dfMeta <- cbind(df_out,subsetMeta)
            
            pcaWithLM <- pca
            summary(pca)
          }
      }
    }
  
    
    nIndividuals <- length(unique(subsetMeta$animal_id))
    nGenes <- dim(filterCountsDF)[1]

    p <- ggplot(dfMeta,aes(x=PC1,y=PC2,color=Trt)) + geom_point(size=3) + 
      xlab(paste0("PC1: ",100 * round(summary(pca)$imp[2,1],3),"% variance")) + 
      ylab(paste0("PC2: ",100 * round(summary(pca)$imp[2,2],3),"% variance")) +
      ggtitle(paste0(gatherSubset,", model:",model), 
              subtitle = paste0(nIndividuals," individuals, ",nGenes," genes",titleText)) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    p
    
    ggsave(plot = p,
           filename = paste0(figDir,gatherSubset,"/pca_LM_",lmUsed,"_model",model,".png"),
           width = 6, height = 4)
    
    p + geom_label_repel(aes(label = subsetMeta$sampleID), 
                         box.padding = 0.35, point.padding = 0.8, 
                         segment.color = 'grey50', force = 10, size = 2)
    
    ggsave(plot = p,
           filename = paste0(figDir,gatherSubset,"/pca_LM_",lmUsed,"_model",model,"_wLabels.png"),
           width = 6, height = 4)
    
    p <- ggplot(dfMeta,aes(x=PC1,y=PC2,color=seqBatch)) + geom_point(size=3) + 
      xlab(paste0("PC1: ",100 * round(summary(pca)$imp[2,1],3),"% variance")) + 
      ylab(paste0("PC2: ",100 * round(summary(pca)$imp[2,2],3),"% variance")) +
      ggtitle(paste0(gatherSubset,", model:",model), 
              subtitle = paste0(nIndividuals," individuals, ",nGenes," genes",titleText)) + 
      theme(plot.title = element_text(hjust = 0.5))
    
    ggsave(plot = p,
           filename = paste0(figDir,gatherSubset,"/pca_LM_",lmUsed,"_model",model,"_colBatch.png"),
           width = 6, height = 4)
    
    ### bring in plots showing secondary correlations
    ## traits is trait matrix of traits of interest (rows are samples, columns are traits)
    
    traitList <- c("Reads","MappedReads","FeatureReads","percentMapped","PercentInFeat","UniqueGenes",
                   "qubitConc","group","libraryBatch","cleanBatch","sampleMonth","climSeason","avgWeight60d","weightSlope","drawTimeHR")
    
    if (length(unique(subsetMeta$seqBatch)) > 1) {
      traitList <- c(traitList,"seqBatch")
    }
    
    if (length(unique(subsetMeta$pregStatusBackDate)) > 1) {
      traitList <- c(traitList,"pregStatusBackDate35")
    }
    
    traits <- subsetMeta[which(colnames(subsetMeta) %in% traitList)]
    
    rownames(traits) <- subsetMeta$sampleID
    
    loads <- pca$x[as.character(rownames(traits)),]
    loads <- cbind(loads,traits)
    
    ## look at effect (p-val from linear model) of each variable on each of top n PCs to 
    nPCs <- min(10,dim(pca$x)[2] - 1)
    pcaDim <- dim(pca$x)[2]
    totalVars <- dim(traits)[2]
    
    varCorrs <- data.frame()
    
    #generates the p-value of the linear model between PC-X and the given trait
    #as an example, PC9 and trait 3
    #-log10(summary(lm(pcLoadings[,9] ~ loads[,pcaDim + 3]))$coef[2,4])
    
    for (v in 1:totalVars) {
      tmpDF <- as.data.frame(cbind(-log10(apply(loads[,1:nPCs],2,function(x){summary(lm(x~loads[,pcaDim + v]))$coef[2,4]})),1:nPCs))
      tmpDF$var <- colnames(loads)[pcaDim + v]
      varCorrs <- rbind(varCorrs, tmpDF)
    }
    
    varCorrs
    
    techCorPlot <- ggplot() +
      geom_hline(yintercept = 2, alpha = .5) +
      geom_point(aes(x = varCorrs$V2, y = varCorrs$V1, col = varCorrs$var), size = 1) +
      xlab("PC Number") + ylab("-log10(p-value)") +
      geom_text_repel(aes(x = subset(varCorrs, V1 > 2)$V2,
                           y = subset(varCorrs, V1 > 2)$V1,
                           label = subset(varCorrs, V1 > 2)$var), 
                       box.padding = 0.35, point.padding = 0.8, 
                       segment.color = 'grey50', force = 10, size = 1.5) +
      labs(col = "Variable") +
      ggtitle(paste0("Variable correlation w PCs"), 
              subtitle = paste0(gatherSubset,", model:",model,titleText)) + 
      theme(plot.title = element_text(hjust = 0.5),
            legend.text = element_text(size=6),
            legend.title = element_text(size = 8),
            legend.key.size = unit(.25, "cm"))
    
    ggsave(plot = techCorPlot,
           filename = paste0(figDir,gatherSubset,"/covariates_LM_",lmUsed,"_model",model,".png"),
           width = 6, height = 4)
    

    techCorPlot <- techCorPlot + ylim(c(0,10))
    
    ggsave(plot = techCorPlot,
           filename = paste0(figDir,gatherSubset,"/covariates_LM_",lmUsed,"_model",model,"_ylim10.png"),
           width = 6, height = 4)

    if (lmMethod == "used") {
      
      #write out counts and residuals to file
      saveRDS(object = residuals, file = paste0(outDir,"/",gatherSubset,"_model",model,"_resids.RDS"))
      saveRDS(object = dfMeta, file = paste0(outDir,"/",gatherSubset,"_model",model,"_dfMeta.RDS"))
      
    }
  }
}




