## ---------------------------
##
## Script name: gatherData.R
##
## Purpose of script: pull together meerkat data for a specific subset of analysis
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

# 
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
#ReSampVar <- "preg"

gatherData <- function(workDir, batchList, focalSex, focalTreatment, minGenes = 10000, 
                       minLogCPM = 5, singleTrt = FALSE, combineTrt = FALSE,
                       longReSample = FALSE, ReSampVar = "status"){
  
  dataDir <- paste0(workDir,"data/")
  outDir <- paste0(workDir,"output/")
  figDir <- paste0(workDir,"figures/")
  
  if (combineTrt) {
    focalTreatment <- "allTrt"
    #override single treatment with allTreatment
    singleTrt <- FALSE
  }
  
  if (singleTrt) {
    gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"only_",minGenes,"genes_",minLogCPM,"lCPM")
    if (longReSample) {
      gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"only_",minGenes,"genes_",minLogCPM,"lCPM_",ReSampVar,"ReSamp")
    }
  } else {
    gatherSubset <- paste0(batchList,"_",focalSex,"_",focalTreatment,"_",minGenes,"genes_",minLogCPM,"lCPM")
  }

  metaBatch1 <- readRDS(paste0(dataDir,"metadata_Meerkat_Batch1.rds"))
  metaBatch2 <- readRDS(paste0(dataDir,"metadata_Meerkat_Batch2.rds"))
  
  ### VVF129 needs to be VVHF129 - was sampled with VVHF127, VVHF132, VVHF135
  
  ### VWM157 hould be VLM157 - according to meerkat database capture 15349 (correct date), the animal
  ### in question is in whiskers, but is actually VLM not VWM
  
  
  #metadata
  #batch 1 fixes
  metaBatch1$Reads <- metaBatch1$inputReads
  metaBatch1$MappedReads <- metaBatch1$mappedReads
  metaBatch1$percentMapped <- metaBatch1$percentMap
  metaBatch1$treatment.x <- metaBatch1$treatment
  
  ### new columns
  metaBatch1$seqBatch <- "UChicB1"
  metaBatch1$Accel <- FALSE
  
  #convert to date
  metaBatch1$sampleDate <- as.Date(metaBatch1$sampleDDMMYY,"%m/%d/%y")
  metaBatch1$sampleYear <- format(metaBatch1$sampleDate, format = "%Y")
  metaBatch1$sampleMonth <- format(metaBatch1$sampleDate, format = "%m")
  
  
  #batch 2 fixes
  metaBatch2$seqBatch <- ifelse(metaBatch2$seqLoc == "UChicago", "UChicB2","Michigan")
  
  metaBatch2$draw_time[which(metaBatch2$draw_time == "'0745")] <- "0745"
  
  jointCols <- intersect(colnames(metaBatch1),colnames(metaBatch2))
  
  countsBatch1 <- readRDS(paste0(dataDir,"rawCounts_Meerkat_Batch1.rds"))
  countsBatch2 <- readRDS(paste0(dataDir,"rawCounts_Meerkat_Batch2.rds"))
  
  if (batchList == "allBatches") {
    #counts
    combinedGeneList <- intersect(rownames(countsBatch1),rownames(countsBatch2))
    countsBatch1 <- countsBatch1[rownames(countsBatch1) %in% combinedGeneList,]
    countsBatch2 <- countsBatch2[rownames(countsBatch2) %in% combinedGeneList,]
    countsDF <- cbind(countsBatch1,countsBatch2)
    
    metaBatch1combine <- metaBatch1[,colnames(metaBatch1) %in% jointCols]
    metaBatch2combine <- metaBatch2[,colnames(metaBatch2) %in% jointCols]
    
    metadata <- rbind(metaBatch1combine,metaBatch2combine)
  } else if (batchList == "Batch2020") {
    combinedGeneList <- rownames(countsBatch1)
    countsDF <- countsBatch1
    
    metaBatch1combine <- metaBatch1[,colnames(metaBatch1) %in% jointCols]
    metadata <- metaBatch1combine
    
  } else if (batchList == "Batch2021") {
    combinedGeneList <- rownames(countsBatch2)
    countsDF <- countsBatch2
    
    metaBatch2combine <- metaBatch2[,colnames(metaBatch2) %in% jointCols]
    metadata <- metaBatch2combine
  } 
  
  #add a few basic bits of info for filtering
  metadata <- metadata[order(metadata$sampleID),]
  countsDF <- countsDF[,order(colnames(countsDF))]
  
  all(metadata$sampleID == colnames(countsDF))
  
  #add percent in feature & unique genes
  metadata$PercentInFeat <- 0
  metadata$ReadsInFeat <- 0
  metadata$UniqueGenes <- 0
  
  for (n in 1:dim(metadata)[1]) {
    metadata$ReadsInFeat[n] <- sum(countsDF[,n])
    metadata$UniqueGenes[n] <- sum(countsDF[,n] != 0)
  }
  
  metadata$PercentInFeat <- metadata$ReadsInFeat / metadata$Reads
  
  #bring in last bit of metadata
  samplingRationale <- read.csv2(paste0(dataDir,"KMPSamplingRationale.csv"), sep = ",")
  birthdays <- read.csv2(paste0(dataDir,"KMPBirthdays.csv"), sep = ",")
  pregLact <- read.csv2(paste0(dataDir,"pregStatus.csv"), sep = ",")
  allPregs <- read.csv2(paste0(dataDir,"KMPPregnancies.csv"), sep = ",")
  allWeights <- read.table(paste0(outDir,"DataSet_Weights_202200714.csv"), sep = ",", header = T)
  
  metadataRationale <- merge(metadata,samplingRationale, by = "sampleID", suffixes = c("",".long"))
  
  ### add birthday in, but not via merge, have to write into each column
  metadataRationale$kmpID <- ""
  metadataRationale$finalBirthday <- ""
  metadataRationale$sampleTypo <- FALSE
  metadataRationale$correctsampleID <- ""
  metadataRationale$finalAgeAtSample <- 0
  
  #all(metadata$animal_id %in% birthdays$sampleID)
  
  for (n in 1:dim(metadataRationale)[1]) {
    bdRow <- which(birthdays$sampleID == metadataRationale$animal_id[n])
    
    metadataRationale$finalBirthday[n] <- birthdays$birthday[bdRow]
    metadataRationale$kmpID[n] <- birthdays$kmpIndividID[bdRow]
    metadataRationale$sampleTypo[n] <- birthdays$sampleTypo[bdRow]
    metadataRationale$correctsampleID[n] <- birthdays$correctSampleID[bdRow]
  }
  
  # sample and first bday is MM / DD / YY
  metadataRationale$birthdateDDMMYY[1]
  metadataRationale$birthdateDDMMYY <- as.Date(metadataRationale$birthdateDDMMYY, format = "%m/%d/%y")
  
  metadataRationale$sampleDDMMYY[1]
  metadataRationale$sampleDDMMYY <- as.Date(metadataRationale$sampleDDMMYY, format = "%m/%d/%y")
  
  #final Birthday is MM / DD / YYYY
  metadataRationale$finalBirthday[1]
  metadataRationale$finalBirthday <- as.Date(metadataRationale$finalBirthday, format = "%m/%d/%y")
  
  #subtract dates to get days, divide by 365 to get years
  metadataRationale$finalAgeAtSample <- as.numeric(metadataRationale$sampleDDMMYY - metadataRationale$finalBirthday) / 365
  metadataRationale$ageAtSample <- as.numeric(metadataRationale$ageAtSample)
  
  metadataRationale$ageAtSample <- metadataRationale$finalAgeAtSample
  metadataRationale$birthdate <- metadataRationale$finalBirthday
  metadataRationale$sampleDate <- metadataRationale$sampleDDMMYY
  
  ### add weights in, but not via merge, have to write into each column
  metadataRationale$avgWeight60d <- 0
  metadataRationale$weightAtSampling <- 0
  metadataRationale$weightChange <- 0
  metadataRationale$weightSlope <- 0
  metadataRationale$weightInt <- 0
  metadataRationale$totalWeights <- 0
  metadataRationale$avgWeight6mo <- 0
  
  #all(metadata$animal_id %in% birthdays$sampleID)
  
  for (n in 1:dim(metadataRationale)[1]) {
    #need both the individual and the sample date
    idRow <- which(allWeights$animal_id == metadataRationale$correctsampleID[n] | allWeights$animal_id == metadataRationale$animal_id[n])
    sampdateRow <- which(allWeights$sampleDate == metadataRationale$sampleDDMMYY[n])
    
    wghtRow <- idRow[which(idRow %in% sampdateRow)]
    
    metadataRationale$avgWeight60d[n] <- allWeights$avgWeight60d[wghtRow]
    metadataRationale$weightAtSampling[n] <- allWeights$weightAtSampling[wghtRow]
    metadataRationale$weightChange[n] <- allWeights$weightChange[wghtRow]
    metadataRationale$weightSlope[n] <- allWeights$weightSlope[wghtRow]
    metadataRationale$weightInt[n] <- allWeights$weightInt[wghtRow]
    metadataRationale$totalWeights[n] <- allWeights$totalWeights[wghtRow]
    metadataRationale$avgWeight6mo[n] <- allWeights$avgWeight6mo[wghtRow]
    
  }
  
  #make draw time a real value
  metadataRationale$drawHour <- floor(as.numeric(metadataRationale$draw_time)/100)
  metadataRationale$drawMinFraction <- (as.numeric(metadataRationale$draw_time) - 100 * metadataRationale$drawHour) / 60
  
  metadataRationale$drawTimeHR <- metadataRationale$drawHour + metadataRationale$drawMinFraction

  
  #remove extraneous columns
  extraCols <- which(colnames(metadataRationale) %in% c("finalAgeAtSample",
                                                        "finalBirthday",
                                                        "date_sampled",
                                                        "sampleDayNum",
                                                        "sampleDDMMYY",
                                                        "sampleDayExcel",
                                                        "birthdateDDMMYY",
                                                        "drawHour",
                                                        "drawMinFraction"))
  
  metadata <- metadataRationale[,-extraCols]
  
  metadata$collectionRationale <- ifelse(metadata$rationale == "U", "population", "longitudinal")
  metadata$qubitConc <- as.numeric(metadata$qubitConc)
  
  ### add pregnancy statuses
  metadata$pregStatus <- FALSE
  metadata$lactStatus <- FALSE
  
  #for every listing in the pregnancy | lactation - flip the appropriate binary variable to TRUE
  #row <- 3
  for (row in 1:dim(pregLact)[1]) {
    pregAnimal <- pregLact$animal_id[row]
    pregDate <- as.Date(pregLact$sampleDate[row], format = "%m/%d/%y")
    
    pregRows <- intersect(which(metadata$animal_id == pregAnimal),
                          which(metadata$sampleDate == pregDate))
    
    if (pregLact$pregLact[row] == "preg") {
      metadata$pregStatus[pregRows] <- TRUE
    } else {
      metadata$lactStatus[pregRows] <- TRUE
    }
  }
  
  ### broader definition of pregnancy
  metadata$pregStatusBackDate <- FALSE
  metadata$pregStatusBackDate35 <- FALSE
  metadata$daysPrePreg <- 0
  metadata$pregResult <- ""
  
  #for each pregnancy, check if a sample in this data falls into that window
  #
  # 44 not pregnant
  # 81 should be
  for (row in 1:dim(allPregs)[1]) {
    pregAnimal <- allPregs$FemaleID[row]
    firstPregDate <- as.Date(allPregs$EarliestPregDate[row], format = "%d-%b-%Y")
    lastPregDate <- as.Date(allPregs$PregEndDate[row], format = "%d-%b-%Y")
    
    firstPregDate35 <- as.Date(allPregs$FirstPregnant[row], format = "%d-%b-%Y") - 35
    
    #which animals were sampled in this window
    animalsSampled <- unique((subset(metadata, sampleDate >= firstPregDate & sampleDate <= lastPregDate)$kmpID))
    
    
    #if the animal was sampled while in the 70 day window:
    if (as.character(pregAnimal) %in% animalsSampled) {
      #find out which metadata rows match this animal and date
      pregRowDate <- intersect(which(metadata$sampleDate >= firstPregDate),
                               which(metadata$sampleDate <= lastPregDate))
      pregRows <- intersect(which(metadata$kmpID == pregAnimal),
                               pregRowDate)
      
      
      #swap those rows to TRUE
      metadata$pregStatusBackDate[pregRows] <- TRUE
      metadata$daysPrePreg[pregRows] <- as.numeric(lastPregDate - metadata$sampleDate[pregRows[1]])
      metadata$pregResult[pregRows] <- allPregs$Outcome[row]
    }
    
    animalsSampled35 <- unique((subset(metadata, sampleDate >= firstPregDate35 & sampleDate <= lastPregDate)$kmpID))
    if (as.character(pregAnimal) %in% animalsSampled) {
      #find out which metadata rows match this animal and date
      pregRowDate <- intersect(which(metadata$sampleDate >= firstPregDate35),
                               which(metadata$sampleDate <= lastPregDate))
      pregRows <- intersect(which(metadata$kmpID == pregAnimal),
                            pregRowDate)
      
      
      #swap those rows to TRUE
      metadata$pregStatusBackDate35[pregRows] <- TRUE
      #metadata$daysPrePreg[pregRows] <- as.numeric(lastPregDate - metadata$sampleDate[pregRows[1]])
      #metadata$pregResult[pregRows] <- allPregs$Outcome[row]
    }
  }
  
  #### FLIP THE PREG BETWEEN 35 DAYS and 70 HERE
  
  #metadata$pregLact <- ifelse(metadata$pregStatusBackDate, "preg",
  #                            ifelse(metadata$lactStatus, "lact","none"))

  metadata$pregLact <- ifelse(metadata$pregStatusBackDate35, "preg",
                              ifelse(metadata$lactStatus, "lact","none"))
  
  #Add weights
  
  ### FILTER DOWN THE DATA TO MATCH THE MODEL SUBSET
  #trim down sampleSubset object with filters and intersect() function
  sampleSubset <- metadata$sampleID
  
  
  #no juvenile meerkats
  adultBinary <- as.numeric(metadata$ageAtSample) > .95 | is.na(as.numeric(metadata$ageAtSample))
  sampleSubset <- intersect(metadata$sampleID[adultBinary],sampleSubset)
  
  #drop Michigan Samples
  sampleSubset <- intersect(metadata$sampleID[metadata$seqBatch != "Michigan"],sampleSubset)
  
  ### set filters
  
  #filter metadata
  sampleSubset <- intersect(subset(metadata, UniqueGenes > minGenes)$sampleID, sampleSubset)
  
  #clean treatment column
  metadata$Trt <- ifelse(metadata$treatment.x == "Control", "Ctrl",
                             ifelse(metadata$treatment.x == "Dex_1uM", "Dex",
                                    ifelse(metadata$treatment.x == "Gard_1ug/mL", "Gard",
                                           ifelse(metadata$treatment.x == "LPS_10ng/mL", "LPS",NA))))
  
  #fix lib and clean batch columns
  metadata$seqLibBatch <- paste0(metadata$seqBatch,"_",metadata$libraryBatch)
  
  #add seasonality from Month
  #https://www.nature.com/articles/nature17986
  metadata$climSeason <- ifelse(metadata$sampleMonth %in% c("10","11","12","01","02","03","04"), "warm-wet",
                                ifelse(metadata$sampleMonth %in% c("05","06","07","08","09"), "cold-dry", "none"))
  
  
  ### only focal treatment
  if (combineTrt) {
    #keep everything
    trtKeep <- unique(metadata$Trt)[!is.na(unique(metadata$Trt))]
    sampleSubset <- intersect(subset(metadata, Trt %in% trtKeep)$sampleID, sampleSubset)
  } else {
    if (singleTrt) {
      sampleSubset <- intersect(subset(metadata, Trt %in% c(focalTreatment))$sampleID, sampleSubset)
    } else {
      sampleSubset <- intersect(subset(metadata, Trt %in% c(focalTreatment,"Ctrl"))$sampleID, sampleSubset)
    }
  }
  
  if (focalSex == "F") {
    focalSexList <- "F"
  } else if (focalSex == "M") {
    focalSexList <- "M" 
  } else {
    focalSexList <- c("M","F")
  }
  
  sampleSubset <- intersect(subset(metadata, sex %in% focalSexList)$sampleID, sampleSubset)
  
  #outlier list - these samples all have a weird logCPM density curve
  outlierSamples <- c("SSUR0335","SSUR0342","SSUR0339","SSUR0340")
  #dex outliers from PCA (high PCA1)
  outlierSamples <- c(outlierSamples,"SSUR2532","SSUR2611","SSUR2614","SSUR2629")
  #gard outlier, same individual
  outlierSamples <- c(outlierSamples,"SSUR2613")
  #single library samples, correlating with PC6
  outlierSamples <- c(outlierSamples,"SSUR0714","SSUR1629")
  #drop outliers
  sampleSubset <- sampleSubset[! sampleSubset %in% outlierSamples]
  
  
  #use the sample IDs to filter the counts object
  countsFiltered <- countsDF[,colnames(countsDF) %in% sampleSubset]
  metaFiltered <- subset(metadata, sampleID %in% sampleSubset)
  
  #arrange both in matching order
  metaFiltered <- metaFiltered[order(metaFiltered$sampleID),]
  countsFiltered <- countsFiltered[,order(colnames(countsFiltered))]
  
  #check our order
  all(metaFiltered$sampleID == colnames(countsFiltered))
  
  ### Center weight, center weight by status
  metaFiltered$weight <- ifelse(metaFiltered$avgWeight60d > 0, metaFiltered$avgWeight60d, metaFiltered$avgWeight6mo)
  #subsetMeta$centWeight
  #subsetMeta$statCentWeight
  
  #unique weights
  weightsOnly <- unique(metaFiltered[,colnames(metaFiltered) %in% c("animal_id","weight","status_dom_sub")])
  
  #center all weights
  weightsOnly$centWeight <- scale(weightsOnly$weight, center = T, scale = F)
  
  #center weight within status
  weightsOnly$statCentWeight <- 0
  weightsOnly$statCentWeight[weightsOnly$status_dom_sub == "S"] <- scale(weightsOnly$weight[weightsOnly$status_dom_sub == "S"], center = T, scale = F)
  weightsOnly$statCentWeight[weightsOnly$status_dom_sub == "D"] <- scale(weightsOnly$weight[weightsOnly$status_dom_sub == "D"], center = T, scale = F)
  
  metaFiltered$centWeight <- 0
  metaFiltered$statCentWeight <- 0
  
  for (r in 1:dim(metaFiltered)[1]) {
    animalID <- metaFiltered$animal_id[r]
    exactWeight <- metaFiltered$weight[r]
    status <- metaFiltered$status_dom_sub[r]
    metaFiltered$centWeight[r] <- as.numeric(subset(weightsOnly, status_dom_sub == status & animal_id == animalID & weight == exactWeight)$centWeight)
    metaFiltered$statCentWeight[r] <- as.numeric(subset(weightsOnly, status_dom_sub == status & animal_id == animalID & weight == exactWeight)$statCentWeight)
  }
  
  subsetMeta <- metaFiltered
  subsetCountsDF <- countsFiltered
  
  allGenes <- rownames(subsetCountsDF)
  
  if (longReSample) {
    if (ReSampVar == "status") {
      mrktCaptures <- unique(subsetMeta[,colnames(subsetMeta) %in% c("correctsampleID","captureref","status_dom_sub")])
      #in here twice
      dupSamples <- names(table(mrktCaptures$correctsampleID)[table(mrktCaptures$correctsampleID) > 2])
      dropCaptures <- vector()
      
      #s <- "VLF230"
      for (s in dupSamples) {
        if ( length(unique(subset(mrktCaptures, correctsampleID == s)$status_dom_sub)) > 1 ) {
          #add the sub samples to a list to drop, by capture ref
          dropCaptures <- c(dropCaptures,
                            subset(mrktCaptures, correctsampleID == s & status_dom_sub == "S")$captureref)
        }
      }
    } else if (ReSampVar == "preg") {
      mrktCaptures <- unique(subsetMeta[,colnames(subsetMeta) %in% c("correctsampleID","captureref","pregStatusBackDate35")])
      #in here twice
      dupSamples <- names(table(mrktCaptures$correctsampleID)[table(mrktCaptures$correctsampleID) > 2])
      dropCaptures <- vector()
      
      #s <- "VLF230"
      for (s in dupSamples) {
        if ( length(unique(subset(mrktCaptures, correctsampleID == s)$pregStatusBackDate35)) > 1 ) {
          #add the sub samples to a list to drop, by capture ref
          dropCaptures <- c(dropCaptures,
                            subset(mrktCaptures, correctsampleID == s & pregStatusBackDate35 == FALSE)$captureref)
        }
      }
    }
    
    subsetMeta <- subset(subsetMeta, ! captureref %in% dropCaptures)
    subsetCountsDF <- subsetCountsDF[,colnames(subsetCountsDF) %in% subsetMeta$sampleID]
  }
  

  #print some facts about the dataset (number of d/s etc) to file
  sink(paste0(outDir,"/",gatherSubset,"_sampleDescription.txt"))

  print(paste0("There are ",length(unique(subsetMeta$animal_id))," unique meerkats sampled"))
  table(unique(subsetMeta[,c(2,5)])[,2])
  table(unique(subsetMeta[,c(2,6)])[,2])
  uniqAnimals <- unique(subsetMeta[,c(which(colnames(subsetMeta) %in% c("animal_id","status_dom_sub","sex")))])
  
  table(paste0(uniqAnimals$sex,"_",uniqAnimals$status_dom_sub))
  print(paste0("With ",length(unique(subsetMeta$captureref))," distinct captures"))
  
  print(paste0("and ",dim(subsetMeta)[1]," total challenge samples"))
  print(paste0("There are ",length(unique(subsetMeta$transitionID)) - 1," unique meerkat transitions"))
  
  print(paste0("With ",dim(subset(subsetMeta, rationale %in% c("newDom","conDom")))[1]," dominant samples"))
  print(paste0("and ",dim(subset(subsetMeta, rationale %in% c("subMale","subFemale","conSub")))[1]," subordinate samples"))

  sink()
  
  #?#?# check the metadata for correlations!
  corCheck <- c("group","status_dom_sub","ageAtSample","FeatureReads","seqLibBatch","qubitConc",
                "Reads","MappedReads","percentMapped","Trt","sampleYear","sampleMonth","PercentInFeat",
                "UniqueGenes","climSeason","lactStatus","pregStatusBackDate","pregLact",
                "avgWeight60d","weightSlope","drawTimeHR")
  
  if (focalSex != "F" & focalSex != "M") {
    corCheck <- c(corCheck,"sex")
  }
  
  if (batchList == "allBatches") {
    corCheck <- c(corCheck,"seqBatch")
  }
  
  metaCorrCheck <- subsetMeta[,colnames(subsetMeta) %in% corCheck]
  
  #FOR EACH COMBO, ggplot them against each other and save
  corrCheckDir <- paste0(figDir,gatherSubset,"/corrCheck/")
  dir.create(paste0(figDir,gatherSubset))
  dir.create(corrCheckDir)
  
  nVars <- dim(metaCorrCheck)[2]
  
  
  for (i in 1:nVars) {
    for (j in 1:nVars) {
      if (i != j) {
        if (j > i) {
          #x is i variable
          xJitter <- .5
          
          if (class(metaCorrCheck[,i]) == "numeric") {
            xJitter <- .05 * (range(metaCorrCheck[!is.na(metaCorrCheck[,i]),i])[2] - 
                                range(metaCorrCheck[!is.na(metaCorrCheck[,i]),i])[1])
          } 
          
          #y is j variable
          yJitter <- .5
          
          if (class(metaCorrCheck[,j]) == "numeric") {
            yJitter <- .05 * (range(metaCorrCheck[!is.na(metaCorrCheck[,j]),j])[2] - 
                                range(metaCorrCheck[!is.na(metaCorrCheck[,j]),j])[1])
          } 
          
          #if there is a simple variable, color it
          if (length(unique(metaCorrCheck[,i])) < 6) {
            if (length(unique(metaCorrCheck[,j])) < 6) {
              corPlot <- ggplot() +
                geom_jitter(aes(x = metaCorrCheck[,i],
                                y = metaCorrCheck[,j],
                                col = metaCorrCheck[,i],
                                shape = metaCorrCheck[,j]),
                            alpha = .5,
                            height = yJitter,
                            width = xJitter) +
                xlab(colnames(metaCorrCheck)[i]) +
                ylab(colnames(metaCorrCheck)[j]) +
                labs(col = colnames(metaCorrCheck)[i],
                     shape = colnames(metaCorrCheck[j]))
            } else {
              corPlot <- ggplot() +
                geom_jitter(aes(x = metaCorrCheck[,i],
                                y = metaCorrCheck[,j],
                                col = metaCorrCheck[,i]),
                            alpha = .5,
                            height = yJitter,
                            width = xJitter) +
                xlab(colnames(metaCorrCheck)[i]) +
                ylab(colnames(metaCorrCheck)[j]) +
                labs(col = colnames(metaCorrCheck)[i])
            }
          } else if (length(unique(metaCorrCheck[,j])) < 6) {
            corPlot <- ggplot() +
              geom_jitter(aes(x = metaCorrCheck[,i],
                              y = metaCorrCheck[,j],
                              col = metaCorrCheck[,j]),
                          alpha = .5,
                          height = yJitter,
                          width = xJitter) +
              xlab(colnames(metaCorrCheck)[i]) +
              ylab(colnames(metaCorrCheck)[j]) +
              labs(col = colnames(metaCorrCheck)[j])
          } else {
            corPlot <- ggplot() +
              geom_jitter(aes(x = metaCorrCheck[,i],
                              y = metaCorrCheck[,j]),
                          alpha = .5,
                          height = yJitter,
                          width = xJitter) +
              xlab(colnames(metaCorrCheck)[i]) +
              ylab(colnames(metaCorrCheck)[j])
            
          }
          
          if ( class(metaCorrCheck[,i]) %in% c("integer","numeric") ) {
            if ( class(metaCorrCheck[,j]) %in% c("integer","numeric") ) {
              corPlot <- corPlot + geom_smooth(aes(x = metaCorrCheck[,i],
                                                   y = metaCorrCheck[,j]),
                                               method = "lm",
                                               formula = "y~x",
                                               se = FALSE)
            }
          }
          
          plotTitle <- paste0(colnames(metaCorrCheck)[j],"_v_",colnames(metaCorrCheck)[i],".png")
          
          ggsave(plot = corPlot,
                 filename = paste0(corrCheckDir,plotTitle), width = 6, height = 4)
        }
      }
    }
  }
  
  #write the counts and metadata out
  saveRDS(object = subsetMeta, 
              file = paste0(outDir,"/",gatherSubset,"_metadata.RDS"))
  saveRDS(object = subsetCountsDF, 
              file = paste0(outDir,"/",gatherSubset,"_counts.RDS"))
  
  
  #go ahead and filter the count object | by sig gene or logCPM
  if (minLogCPM == "status") {
    nSamples <- dim(subsetMeta)[1]
    nIndividuals <- length(unique(subsetMeta$animal_id))
    
    sigGenes <- readRDS("~/Dropbox (Personal)/meerkats/RNA/fullAnalysis/data/sigStatusGenes.RDS")
    passLCPM <- sigGenes
    
    passLCPMuniq <- unique(passLCPM)
    nGenes <- length(passLCPMuniq)
    
    modelCountsLCPMfilt <- subsetCountsDF[rownames(subsetCountsDF) %in% passLCPMuniq,]
    
  } else {
    lCPMthresh <- as.numeric(minLogCPM)
    allGenes <- rownames(subsetCountsDF)
    
    nSamples <- dim(subsetMeta)[1]
    nIndividuals <- length(unique(subsetMeta$animal_id))
    
    ## calc avg log CPM by condition (any number of conditions)
    conditionList <- unique(subsetMeta$Trt)
    
    passLCPM <- vector()
    
    for (c in conditionList) {
      geneLCPMcond <- aveLogCPM(subsetCountsDF[,subsetMeta$Trt == c])
      passLCPM <- c(passLCPM,allGenes[geneLCPMcond > lCPMthresh])
    }
    
    passLCPMuniq <- unique(passLCPM)
    nGenes <- length(passLCPMuniq)
    
    ## use list of genes that pass to cut down list
    modelCountsLCPMfilt <- subsetCountsDF[rownames(subsetCountsDF) %in% passLCPMuniq,]
  }
  
  logCPM <- cpm(modelCountsLCPMfilt, log=TRUE)
  plotLCPMfilt <- melt(logCPM)
  ggplot(plotLCPMfilt) +
    geom_density(aes(x = value, col = Var2)) + 
    ggtitle(paste0("Filtered, logCPM values, ",nIndividuals," Individuals, ",nGenes," genes"),
            subtitle = gatherSubset) + 
    theme(legend.position = "none")
  ggsave(paste0(figDir,gatherSubset,"/logCPM_density_plot.png"),
         width = 6, height = 4)
  
  saveRDS(object = modelCountsLCPMfilt, 
          file = paste0(outDir,"/",gatherSubset,"_filteredCounts.RDS"))
  
}
