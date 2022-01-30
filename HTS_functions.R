library(edgeR)

prepare_matrix <- function(file,...){

	#For voom we must have a numeric matrix with
	#row names = gene ids
	#receives counts file, with first column being "entrezIDs"
	#returns matrix for voom

	counts <- read.table(file,stringsAsFactors = F,...)
	#remove duplicate rows (some gene identifiers can map to two entrez ids)
	#get row with max mean
	counts2 <- counts[, 2:ncol(counts)]
	counts2 <- getMaxMean(counts)
	
	#change row names, remove gene id column, turn to matrix
	row.names(counts2) <- counts2$entrezIDs
	counts2 <- counts2[ , colnames(counts2) != "entrezIDs"]
	counts2 <- as.matrix(counts2)
	return(counts2)
}
proccessHTS <- function(counts, designMatrix, log = ""){
  
  # preprocess counts file before feeding to limma using voom
  # receives count matrix and design matrix and outputs
  # and recieves an EList object with normalized counts
  # and weighst of samples inversely proprtional to estimated variance
  # based on counts as per voom methodology
  # this EList object can be inserted as annotatedGSE in generalizedTtest function
  
  #normalization
  d0 <- DGEList(counts)
  d0 <- calcNormFactors(d0)
  cutoff <- 1
  drop <- which(apply(cpm(d0), 1, max) < cutoff)
  d <- d0[-drop,] 

  #reorder the design matrix so row ordering matches column ordering of annotatedGseset
  #limma will not look for rownames
  designMatrix <- designMatrix[order(match(rownames(designMatrix), colnames(counts))), , drop = F]
  
  #weight calculation
  y <- voom(d, designMatrix, plot = F)

  #generalizedTtest is meant to work with arrays too
  #as such expexts a gene field where probes (rows in the fit)
  #are annotated with entrez gene ids, so we need to tell it
  #that our rows are gene ids...which they are
  #wastes a bit of space, but allows to reuse generalizedTtest
  y$genes <- data.frame(entrezIDs = row.names(d))
  row.names(y$genes) <- y$genes$entrezIDs
  cat("VOOM weights calculated, passing to limma\n",file = log, append = T)
  return(y)}
  '
  fit <- lmFit(y, designMatrix)
  fit<-contrasts.fit(fit, contMatrix)
  fit<-eBayes(fit)
  
  for(coeff in colnames(contMatrix)){
    
    tt <-topTable(fit, coef = coeff, number= nrow(y), adjust.method="BH", sort.by="logFC")
    tt$entrezIDs <- rownames(y$E)
    signGenes <- tt[ tt$adj.P.Val < 0.05, c("entrezIDs","logFC"), drop = F]
    write.table(tt, file = paste( direc,"/topTable/",file, "_top_table_", coeff ,sep = "")  )
    
    if(nrow(signGenes) == 0){
      cat(paste("No significant genes after p value adjustment found for", coeff), sep = "\n"
          , file = log, append = T)
      
    }else{
      getBINGSpaceGenes(signGenes = signGenes, probeToEntrezIDs = F, name = paste(file,coeff, sep = "_"),log = log, direc = direc, col = "ID" )
    }
  }
}'
