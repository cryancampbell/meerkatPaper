## 
## ---------------------------
##
## Script name: compLMMBeta.R
##
## Purpose of script: take two emmreml outputs and compare the betas of a given variable
##
## Author: CRC
##
## Date Created: 2022-02-22
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------

 emremlX <- "~/Dropbox (Personal)/meerkats/RNA/fullAnalysis/output/Batch2021_F_LPS_10000genes_4lCPM_modelF_nested_emmreml"
 colnameX <- "treatment"
 Xtitle <- "Batch2021"

 emremlY <- "~/Dropbox (Personal)/meerkats/RNA/fullAnalysis/output/Batch2020_F_LPS_10000genes_4lCPM_modelF_nested_emmreml"
 colnameY <- "treatment"
 Ytitle <- "Batch2020"

 figBetaDir <- paste0(figDir,"betaComp")

 rnaCounts <- readRDS("~/Dropbox (Personal)/meerkats/RNA/fullAnalysis/output/allBatches_F_LPS_10000genes_4lCPM_counts.RDS")
 meanGE <- rowMeans(rnaCounts)
 topNpercent <- 25

 varA <- varB <- "Gard"
 modelA <- modelB <- "modelG"
#
emremlX = paste0("~/Dropbox (Personal)/meerkats/RNA/fullAnalysis/output/allBatches_F_",varA,"only_10000genes_5lCPM_",modelA,"_basicPreg_emmreml")
colnameX = "status"
Xtitle = paste0(varA,"only")
emremlY = paste0("~/Dropbox (Personal)/meerkats/RNA/fullAnalysis/output/allBatches_F_",varB,"only_10000genes_5lCPM_",modelB,"_basicPreg_emmreml")
colnameY = "preg"
Ytitle = paste0(varB,"only")
figBetaDir = paste0(figDir,"betaComp")
meanGE = rowMeans(rnaCounts)
topNpercent = 25
labelOutlier <- FALSE
compareSigGenes <- TRUE

sigThresh <- 0.10
mashComp <- FALSE
useQval <- TRUE

compLMMBeta <- function(figBetaDir, emremlX, colnameX="treatment", Xtitle,
                        emremlY, colnameY="treatment", Ytitle,
                        meanGE, topNpercent, labelOutlier = FALSE, 
                        compareSigGenes = TRUE, sigThresh = 0.10, mashComp = FALSE, useQval = TRUE){
  
  if ( class(emremlX) == "character" ) {
    emremlX <- read.table(emremlX)
  }
  
  if ( class(emremlY) == "character" ) {
    emremlY <- read.table(emremlY)
  }
  
  mashSave <- ""
  if (mashComp == TRUE) {
    mashSave <- "_mash"
  }
  
  sharedGenes <- intersect(rownames(emremlX),rownames(emremlY))
  
  emremlX$geneName <- rownames(emremlX)
  emremlY$geneName <- rownames(emremlY)
  
  compEmreml <- merge(emremlX,emremlY, by = "geneName", suffixes = c(".x",".y"))
  
  if (paste0("beta_",colnameX) %in% colnames(emremlY)) {
    colX <- which(colnames(compEmreml) == paste0("beta_",colnameX,".x"))
    pvalX <- which(colnames(compEmreml) == paste0("pval_",colnameX,".x"))
    varX <- which(colnames(compEmreml) == paste0("var_beta_",colnameX,".x"))
  } else {
    colX <- which(colnames(compEmreml) == paste0("beta_",colnameX))
    pvalX <- which(colnames(compEmreml) == paste0("pval_",colnameX))
    varX <- which(colnames(compEmreml) == paste0("var_beta_",colnameX))
  }

  if (paste0("beta_",colnameY) %in% colnames(emremlX)) {
    colY <- which(colnames(compEmreml) == paste0("beta_",colnameY,".y"))
    pvalY <- which(colnames(compEmreml) == paste0("pval_",colnameY,".y"))
    varY <- which(colnames(compEmreml) == paste0("var_beta_",colnameY,".y"))
  } else {
    colY <- which(colnames(compEmreml) == paste0("beta_",colnameY))
    pvalY <- which(colnames(compEmreml) == paste0("pval_",colnameY))
    varY <- which(colnames(compEmreml) == paste0("var_beta_",colnameY))
  }
  
  
  compBetaPlot <- ggplot(compEmreml, aes_string(x=paste0(colnames(compEmreml)[colX]), 
                                                y=paste0(colnames(compEmreml)[colY]))) +
    geom_point(alpha = .2) +
    geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0(Xtitle," ",colnameX," betas")) +
    ylab(paste0(Ytitle," ",colnameY," betas")) +
    ggtitle(paste0(Ytitle," v. ",Xtitle," ",colnameX," betas"), 
            subtitle = paste0("R=",signif(cor.test(y = compEmreml[,colY], 
                                                   x = compEmreml[,colX])$estimate, digits = 3),
                              " p=",signif(cor.test(y = compEmreml[,colY], 
                                                    x = compEmreml[,colX])$p.value, digits = 3),
                              " nGenes=",dim(compEmreml)[1]))
  ggsave(plot = compBetaPlot, height = 4, width = 6,
         filename = paste0(figBetaDir,"/",Ytitle,colnameY,"_v_",Xtitle,colnameX,mashSave,".png"))
  
  
  geStrengthShared <- meanGE[names(meanGE) %in% sharedGenes]
  if (all(names(geStrengthShared) == compEmreml$geneName)) {
    compEmreml$GE <- geStrengthShared
  } else {
    #if the lists aren't the same, trim both by common genes
    sharedGEGenes <- intersect(names(geStrengthShared),compEmreml$geneName)
    
    geStrengthShared <- geStrengthShared[names(geStrengthShared) %in% sharedGEGenes]
    compEmreml <- compEmreml[rownames(compEmreml) %in% sharedGEGenes,]
    
    compEmreml$GE <- geStrengthShared
  }
  
  
  #top percent cutoff
  geCutoff <- quantile(compEmreml$GE, probs = (100 - topNpercent)/100)
  compEmremlHighGE <- subset(compEmreml, GE > geCutoff)
  
  compBetaPlotHighGE <- ggplot(compEmremlHighGE, aes_string(x=paste0(colnames(compEmreml)[colX]), 
                                                            y=paste0(colnames(compEmreml)[colY]))) +
    geom_point(alpha = .2) +
    geom_smooth(method = "lm", formula = "y~x") +
    xlab(paste0(Xtitle," ",colnameX," betas")) +
    ylab(paste0(Ytitle," ",colnameY," betas")) +
    ggtitle(paste0(Ytitle," v. ",Xtitle," ",colnameX," betas | top ",topNpercent,"% GE"), 
            subtitle = paste0("R=",signif(cor.test(y = compEmremlHighGE[,colY], 
                                                   x = compEmremlHighGE[,colX])$estimate, digits = 3),
                              " p=",signif(cor.test(y = compEmremlHighGE[,colY], 
                                                    x = compEmremlHighGE[,colX])$p.value, digits = 3),
                              " nGenes=",dim(compEmremlHighGE)[1]))
  ggsave(plot = compBetaPlotHighGE, height = 4, width = 6,
         filename = paste0(figBetaDir,"/",Ytitle,colnameY,"_v_",Xtitle,colnameX,"_top",topNpercent,"ge",mashSave,".png"))
  
  if (labelOutlier == TRUE) {
    
    outliers <- compEmremlHighGE[compEmremlHighGE[,colX] < -.25 & compEmremlHighGE[,colY] > .1,]
    
    compBetaPlotHighGElabel <- ggplot() +
      geom_point(aes(x=compEmremlHighGE[,colX],
                     y=compEmremlHighGE[,colY]), alpha = .2) +
      geom_label_repel(aes(x = outliers[,colX],
                           y = outliers[,colY],
                           label = outliers$geneName)) +
      geom_smooth(aes(x=compEmremlHighGE[,colX],
                      y=compEmremlHighGE[,colY]), method = "lm", formula = "y~x") +
      xlab(paste0(Xtitle," ",colnameX," betas")) +
      ylab(paste0(Ytitle," ",colnameY," betas")) +
      ggtitle(paste0(Ytitle," v. ",Xtitle," ",colnameX," betas | top ",topNpercent,"% GE"), 
              subtitle = paste0("R=",signif(cor.test(y = compEmremlHighGE[,colY], 
                                                     x = compEmremlHighGE[,colX])$estimate, digits = 3),
                                " p=",signif(cor.test(y = compEmremlHighGE[,colY], 
                                                      x = compEmremlHighGE[,colX])$p.value, digits = 3),
                                " nGenes=",dim(compEmremlHighGE)[1]))
    
    ggsave(plot = compBetaPlotHighGElabel, height = 4, width = 6,
           filename = paste0(figBetaDir,"/",Ytitle,colnameY,"_v_",Xtitle,colnameX,"_top",topNpercent,"ge_labels",mashSave,".png"))
  }
  
  if (compareSigGenes == TRUE) {
    #use qvalue if TRUE
    if (useQval) {
      compEmreml$sigInX <- qvalue(compEmreml[,pvalX])$qvalue < sigThresh
      compEmreml$sigInY <- qvalue(compEmreml[,pvalY])$qvalue < sigThresh
    } else {
      #create a column for "significant in varX"
      compEmreml$sigInX <- compEmreml[,pvalX] < sigThresh
      #create a column for "significant in varY"
      compEmreml$sigInY <- compEmreml[,pvalY] < sigThresh
    }
    
    #create a column that combines the two "sig in none, X, Y, both"
    compEmreml$sigInBoth <- ifelse(compEmreml$sigInX & compEmreml$sigInY, "both",
                                   ifelse(compEmreml$sigInX, "sigX",
                                          ifelse(compEmreml$sigInY, "sigY", "none")))
    
    compEmreml$sigInBoth <- factor(compEmreml$sigInBoth, levels = c("none","sigX","sigY","both"))
    
    fishDat <- data.frame(
      "sigX" = c(sum(compEmreml$sigInBoth == "both"),sum(compEmreml$sigInBoth == "sigX")),
      "nsX" = c(sum(compEmreml$sigInBoth == "sigY"),sum(compEmreml$sigInBoth == "none")),
      stringsAsFactors = FALSE
    )
    colnames(fishDat) <- c("sigX", "nsX")
    rownames(fishDat) <- c("sigY", "nsY")
    
    
    fishP <- signif(fisher.test(fishDat)$p.value, digits = 3)
    sigGeneDiff <- signif(fishDat - chisq.test(fishDat)$expected, digits = 3)
    excessBoth <- sigGeneDiff[1,1]
    
    compEmremlSort <- compEmreml[order(compEmreml$sigInBoth),]
    
    compBetaPlot <- ggplot(compEmremlSort, aes_string(x=paste0(colnames(compEmremlSort)[colX]), 
                                                  y=paste0(colnames(compEmremlSort)[colY]),
                                                  col = "sigInBoth")) +
      geom_point(alpha = 3/4) +
      #geom_smooth(method = "lm", formula = "y~x") +
      xlab(paste0(Xtitle," ",colnameX," betas")) +
      ylab(paste0(Ytitle," ",colnameY," betas")) +
      ggtitle(paste0(Ytitle," v. ",Xtitle," ",colnameX," betas"), 
              subtitle = paste0("R=",signif(cor.test(y = compEmremlSort[,colY], 
                                                     x = compEmremlSort[,colX])$estimate, digits = 3),
                                " p=",signif(cor.test(y = compEmremlSort[,colY], 
                                                      x = compEmremlSort[,colX])$p.value, digits = 3),
                                " nGenes=",dim(compEmremlSort)[1],
                                "\nFETp=",fishP,
                                " | both obs=",fishDat[1,1],
                                ", exp=",signif(chisq.test(fishDat)$expected[1,1], digits = 3))) +
      scale_color_manual(values = c("gray70",gg_color_hue(4)[c(1,3,4)]))
    
    ggsave(plot = compBetaPlot, height = 4, width = 6,
           filename = paste0(figBetaDir,"/",Ytitle,colnameY,"_v_",Xtitle,colnameX,"_colSigGenes",mashSave,".png"))
    
    if (mashComp == FALSE) {
      compEmreml$stdBetaX <- compEmreml[,colX] / sqrt(compEmreml[,varX])
      compEmreml$stdBetaY <- compEmreml[,colY] / sqrt(compEmreml[,varY])
    
      compBetaPlot <- ggplot(compEmreml, aes_string(x="stdBetaX", 
                                                        y="stdBetaY",
                                                        col = "sigInBoth")) +
        geom_point(alpha = 3/4) +
        #geom_smooth(method = "lm", formula = "y~x") +
        xlab(paste0(Xtitle," ",colnameX," std. betas")) +
        ylab(paste0(Ytitle," ",colnameY," std. betas")) +
        ggtitle(paste0(Ytitle," v. ",Xtitle," ",colnameX," std. betas"), 
                subtitle = paste0("R=",signif(cor.test(y = compEmremlSort[,colY], 
                                                       x = compEmremlSort[,colX])$estimate, digits = 3),
                                  " p=",signif(cor.test(y = compEmremlSort[,colY], 
                                                        x = compEmremlSort[,colX])$p.value, digits = 3),
                                  " nGenes=",dim(compEmremlSort)[1],
                                  "\nFETp=",fishP,
                                  " | both obs=",fishDat[1,1],
                                  ", exp=",signif(chisq.test(fishDat)$expected[1,1], digits = 3))) +
        scale_color_manual(values = c("gray70",gg_color_hue(4)[c(1,3,4)]))
      
      ggsave(plot = compBetaPlot, height = 4, width = 6,
             filename = paste0(figBetaDir,"/",Ytitle,colnameY,"_v_",Xtitle,colnameX,"_stdBeta_colSigGenes.png"))
      
      compBetaPlot <- ggplot(compEmreml, aes_string(x="stdBetaX", 
                                                    y="stdBetaY",
                                                    col = "sigInBoth")) +
        geom_point(alpha = 2/3) +
        #geom_smooth(method = "lm", formula = "y~x") +
        xlab(paste0(Xtitle," ",colnameX," std. betas")) +
        ylab(paste0(Ytitle," ",colnameY," std. betas")) +
        ggtitle(paste0(Ytitle," v. ",Xtitle," ",colnameX," std. betas"), 
                subtitle = paste0("R=",signif(cor.test(y = compEmremlSort[,colY], 
                                                       x = compEmremlSort[,colX])$estimate, digits = 3),
                                  " p=",signif(cor.test(y = compEmremlSort[,colY], 
                                                        x = compEmremlSort[,colX])$p.value, digits = 3),
                                  " nGenes=",dim(compEmremlSort)[1],
                                  "\nFETp=",fishP,
                                  " | both obs=",fishDat[1,1],
                                  ", exp=",signif(chisq.test(fishDat)$expected[1,1], digits = 3))) +
        scale_color_manual(values = c("grey70","red","blue","purple")) +
        theme_minimal() +
        theme(axis.text = element_text(size = 15))
      
      compBetaPlot
      ggsave(plot = compBetaPlot, height = 4, width = 6,
             filename = paste0(figBetaDir,"/",Ytitle,colnameY,"_v_",Xtitle,colnameX,"_stdBeta_colSigGenes_talkFormat.png"))
      
    }
  }
}
