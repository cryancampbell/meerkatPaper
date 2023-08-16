## 
## ---------------------------
##
## Script name: gsea.R
##
## Purpose of script: run GSEA
##
## Author: CRC | via NSM, JAA, & TV
##
## Date Created: 2022-07-15
##
##
## ---------------------------
##
## Notes:
##   
##
## ---------------------------


GSEA=function(effect_size_matrix,gene_set_list,gene_GO_associations, weight=1,num.cores=NA,pvals=F){
  ## effect_size_matrix == gene x permutations matrix where rownames are the genes and the first column is the observed effect sizes (i.e., columns 2:ncol are the permuted betas)
  ##effect_size_matrix = effect_size_matrix; gene_set_list = unique(GO_res$GO);gene_GO_associations = gene.set2;weight = 1; pvals=F
  if (is.na(num.cores)) num.cores=detectCores()
  clus = makeCluster(getOption("cl.cores", num.cores))
  clusterExport(clus,varlist=ls(),envir=environment())
  a=parApply(clus,effect_size_matrix,2,function(tmp_data){
    names(tmp_data)=rownames(effect_size_matrix)
    tmp_data=sort(tmp_data,decreasing = T)
    correl.vector<-tmp_data
    tmp_all_genes<-names(tmp_data)
    #So that we ultimately weight on our betas when weight==1
    weighted.score.type=weight
    Reduce(rbind,lapply(gene_set_list,function(GOcat){
      tmp_set <- gene_GO_associations$Associated.Gene.Name[gene_GO_associations$GO.Term.Accession == GOcat]
      #First we're indexing what genes are or aren't in our gene set (1 vs 0)
      tag.indicator <- sign(match(tmp_all_genes, tmp_set, nomatch=0)) 
      no.tag.indicator <- 1 - tag.indicator 
      N <- length(tmp_all_genes) 
      Nh <- length(tmp_set) 
      Nm <-  N - Nh 
      #This is just if we want to do a KS test(maybe if we want to compare later to some of your results using GAGE); weight==0
      if (weight == 0) {
        correl.vector <- rep(1, N)
      }
      alpha <- weighted.score.type
      if (pvals==T) {correl.vector=-log10(correl.vector)}
      correl.vector <- abs(correl.vector**alpha)
      #Normalizing our Incremental increase by the sum of the abs our "correlations" (sum.correl.tag) and 
      #normalizing our incremental decrease by the number of genes not in our set (Nm, which remember came from N-Nh above).
      sum.correl.tag  <- sum(correl.vector[which(tag.indicator == 1)])
      #attempting to remove the normalization
      norm.tag    <- 1.0 / sum.correl.tag
      norm.no.tag <- 1.0 / Nm
      #Then we just walk down the gene list adding or substracting as such
      ## don't use normalization in the math
      RES <- cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag)
      #RES <- cumsum(tag.indicator * correl.vector - no.tag.indicator)
      #Then obviously this is for maximum positively vs negatively enriched 
      max.ES <- max(RES)
      min.ES <- min(RES)
      if (max.ES > - min.ES) {
        #      ES <- max.ES
        ES <- signif(max.ES, digits = 5)
        arg.ES <- which.max(RES)
      } else {
        #      ES <- min.ES
        ES <- signif(min.ES, digits=5)
        arg.ES <- which.min(RES)
      }
      #if (length(gene_set)==1) {plot(cumsum(tag.indicator * correl.vector * norm.tag - no.tag.indicator * norm.no.tag),type = "l",xlab="Gene Index",ylab="Running Enrichment Score",xaxt="n"); axis(1,at=tag.indicator*1:length(tag.indicator),labels=T,tick=T)}
      return(c(max.ES,min.ES))}))
  })
  ## Returns two columns, one for upreg and one for downreg
  a=data.frame(a)
  a_upreg=a[1:length(gene_set_list),]
  a_downreg=a[(length(gene_set_list)+1):nrow(a),]
  rownames(a_upreg)=names(gene_set_list)
  rownames(a_downreg)=names(gene_set_list)
  #rownames(a)=names(gene.set)
  ## This doesn't seem right. Why would the pvals be lower if the max doesn't exceed the null?
  pvals_upreg=apply(a_upreg,1,function(x){
    #signif((sum(x[1] >= x[2:length(x)])+1)/(length(x)), digits=5)
    ee <- ecdf(x[2:length(x)]); 1 - ee(x[1])
  })
  names(pvals_upreg)=names(gene_set_list)
  pvals_downreg=apply(a_downreg,1,function(x){
    #signif((sum(x[1] <= x[2:length(x)])+1)/(length(x)), digits=5)
    ee <- ecdf(x[2:length(x)]); ee(x[1])
  })
  names(pvals_downreg)=names(gene_set_list)
  pvals=cbind(pvals_upreg,pvals_downreg); colnames(pvals)=c("pval_upreg","pval_downreg")
  output=list(a_upreg,a_downreg,pvals)
}