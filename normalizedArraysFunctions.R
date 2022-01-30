suppressMessages(library(GEOquery))
suppressMessages(library('org.Hs.eg.db'))
suppressMessages(library(limma))
suppressMessages(library(data.table))

BINGspace <- read.table("/home/mgarranzo/linear_models/BINGspace.tsv", header = T, sep = "\t")
get_file_name <- function(filepath){
  #get base filename without extensions
  return(sub(pattern = "(\\..*)+$", replacement = "", basename(filepath))) 
}


annotateGSE <- function(GSEfile, GPLdir = "/tmp", log = "", GPLFile = "", database = "org.Hs.egREFSEQ2EG"){
  
  cat(paste("#",GSEfile, sep= ""), sep = "\n", log = log, append = T)
  eset <- getGEO(filename = GSEfile, AnnotGPL = FALSE, getGPL = FALSE)
  print(length(GPLFile))
  if(GPLFile != ""){
  
	#Column one: Ids
	#column 2: refseq ids
	#header

	mapping <- read.table(GPLFile, sep = "\t", header = T, comment.char = "#")
    	#genebank <- as.vector(gpleset@dataTable@table[[isGeneBank]])
    	entrezIDs <- getEntrezIDS(mapping[ ,2], as.list(eval(parse(text = database))) )
    	annot <- data.frame("entrezIDs" = entrezIDs)
    	rownames(annot) <- mapping[ ,1]
	#print(head(annot))
  }else{
  #print("....")
  	annot <- checkAnnotations(eset@annotation, eset, GPLdir, log = log)
  }
  if(!is.null(annot)){
    eset@featureData@data <- annot
    eset@featureData@data$entrezIDs <- as.character(eset@featureData@data$entrezIDs)
    order_gse <- match(rownames(eset@assayData$exprs), rownames(eset@featureData@data)) #toptable requires same ordering in featureData as in rows...
    eset@featureData@data <- eset@featureData@data[order_gse,,drop = F]

  }

  return(eset)
  
}
checkAnnotations <- function(GPL, eset, GPLdir = "/tmp", log = ""){
  
  #In microarray experiments, the data matrix returns GPL in annotation
  #need to map probeset IDs to gene names
  #check first if there a bioconductor package
  #if not, annotation using GPL info. 
  
  #if (! "GEOmetadb" %in% .packages()){  library("GEOmetadb")}

  #if(!file.exists('/home/marco/Marco_TFM/Marco_TFM/Rscripts/GEOmetadb.sqlite')) {getSQLiteFile()}
  #con <- dbConnect(SQLite(),'/home/marco/Marco_TFM/Marco_TFM/Rscripts/GEOmetadb.sqlite')
  #query <- paste("select bioc_package from gpl where gpl = '", GPL, "'" , sep = "")
  
  #isPackage = dbGetQuery(con, query )[, "bioc_package" ]
  #dbDisconnect(conn = con)
  isPackage <- NA
  if(is.na(isPackage)){
    
    cat(paste("No package found for GPL", GPL, sep = " "), sep = "\n", file = log, append = T)
    probesetIDstoEntrez <- (annotationGPL(GPL, gset, GPLdir, log = log)) #first column, probe IDs, second column, gene ids
    return(probesetIDstoEntrez)
    
  }

  cat(paste("Package found for GPL", GPL, ". Checking for", isPackage, "in bioconductor", sep = " "), sep = "\n", file = log, append = T)
  annotation <- paste(isPackage, ".db", sep = "")
  return(annotationBiocPackages(eset, annotation))
} 



annotationBiocPackages <- function(eset, annotation){
  
  if (! "annotate" %in% .packages()){  library("annotate")}
  if (!annotation %in% rownames(installed.packages())){
    
    BiocManager::install(annotation)
    library(annotation, character.only = TRUE)
  }else if (!annotation %in% .packages()){
    library(annotation, character.only = TRUE)
  }
  GeneSymbols <- getSYMBOL(rownames(eset), annotation)
  entrezIDs <-getEntrezIDS(toupper(GeneSymbols), as.list(org.Hs.egALIAS2EG))
  df <- data.frame("entrezIDs" = entrezIDs, row.names = rownames(eset))

  return(df)
}



annotationGPL <- function(GPL, eset, GPLdir = "/tmp", log = ""){
  
  #First check if available in user specified  directory (default tmp)
  #or else download it using getGEO
  #checking tmp is default behaviour, but if you have a specific directory for this you can use that
  tmpFiles <- list.files(path = GPLdir, recursive = TRUE, full.names = TRUE)
  isFile <- grep(paste(GPL,"(_family)?(.annot|.soft)", sep = ""), tmpFiles, value = TRUE)
  cat(paste("Checking for GPL soft file in", GPLdir, sep = " "), sep = "\n", file = log, append = T)

  if (length(isFile) == 0){

    cat("File not found. Downloading file", sep = "\n", file = log, append = T)
    gpleset <- getGEO(GEO = GPL, destdir = GPLdir)
  }
  else{
    
    cat(paste("Using", isFile[1], sep = " "), sep = "\n", file = log, append = T)
    gpleset <- getGEO(filename = isFile[1], AnnotGPL = FALSE, getGPL = FALSE)
  }
  
  #just in case, some GPLS have an NA....
  theNa <- is.na(gpleset@dataTable@table$ID)
  nNa <- sum(theNa)
  if(nNa > 0){
    gpleset@dataTable@table[theNa, ]$ID <- paste("NA", 1:sum(theNa), sep = "") 
  }

  #Here comes the tricky part:
  #GPL fields per file can vary. Here are some rules for looking for
  #fields which might be of interest for looking for gene symbols
  #We look for entrez, genbank or ensembl, 
  #since they tend to be the most "polished" fields, if present
  #search is not  exhaustive, just ad-hoc rules that may or may not work
  #if it doesnt, look at the file manually and write a custom script
  #additionally, more than one field can have a description which matches
  #our interest. In that case we always choose the first time
  
  #first just look for genebank ids
  isEntrez <- grep("(entrez( |_))?gene( |_)id" , gpleset@dataTable@columns$Description, ignore.case=TRUE)
  isUni <- grep("uni",gpleset@dataTable@columns$Description, ignore.case = TRUE)
  isEnsemBl <- grep("ensembl", gpleset@dataTable@columns$Description, ignore.case=TRUE)
  isEntrez <- setdiff(isEntrez,isUni) 
  isEntrez <- setdiff(isEntrez,isEnsemBl) 
  
  if(length(isEntrez) != 0){
 
    isEntrez <- isEntrez[1]
    cat(paste("Entrez identifiers found in column", isEntrez, ":",
                colnames(gpleset@dataTable@table)[isEntrez], "(",
                gpleset@dataTable@columns$Description[isEntrez] ,")"), sep = "\n"
                , file = log, append = T)
    
    entrezIDs <- gpleset@dataTable@table[isEntrez]
    rownames(entrezIDs) <- gpleset@dataTable@table$ID
    colnames(entrezIDs) <- "entrezIDs"
    return(entrezIDs)
  }
  
  isRefseq <- grep("(refseq)" , gpleset@dataTable@columns$Description, ignore.case=TRUE)
  if(length(isRefseq) != 0){
    
    isRefseq <-  isRefseq[1]  
    cat(paste("RefSeq identifiers found in column", isRefseq, ":",
                colnames(gpleset@dataTable@table)[isRefseq], "(",
          gpleset@dataTable@columns$Description[isRefseq] ,")"), sep = "\n"
          , file = log, append = T)
    
    refseq <- as.vector(gpleset@dataTable@table[[isRefseq]])
    entrezIDs <- getEntrezIDS(refseq, as.list(org.Hs.egREFSEQ2EG))
    df <- data.frame("entrezIDs" = entrezIDs)
    rownames(df) <- gpleset@dataTable@table$ID
    return(df)
  }
  

  #One field can have one name in one GPL and another one in another GPL
  #can map entrez gene ids to many other things, and viceversa
  isGeneBank <- grep("GB_.*" ,gpleset@dataTable@columns$Column)
  if(length(isGeneBank) != 0){
    
    isGeneBank <- isGeneBank[1] 
    cat(paste("Gene bank identifiers found in column", isGeneBank, ":",
                colnames(gpleset@dataTable@table)[isGeneBank], "(",
          gpleset@dataTable@columns$Description[isGeneBank] ,")"), sep = "\n"
          , file = log, append = T)

    genebank <- gpleset@dataTable@table[[isGeneBank]]
    entrezIDs <- getEntrezIDS(genebank, as.list(org.Hs.egACCNUM2EG))
    df <- data.frame("entrezIDs" = entrezIDs)
    rownames(df) <- gpleset@dataTable@table$ID
    return(df)
    
  }

  
  isGeneBank <- grep("(genbank|gene bank|genebank)" , gpleset@dataTable@columns$Description, ignore.case=TRUE)
  if(length(isGeneBank) != 0){
    
    isGeneBank <-  isGeneBank[1]  
    cat(paste("Gene bank identifiers found in column", isGeneBank, ":",
                colnames(gpleset@dataTable@table)[isGeneBank], "(",
          gpleset@dataTable@columns$Description[isGeneBank] ,")"), sep = "\n"
          , file = log, append = T)
    
    genebank <- as.vector(gpleset@dataTable@table[[isGeneBank]])
    entrezIDs <- getEntrezIDS(genebank, as.list(org.Hs.egACCNUM2EG))
    df <- data.frame("entrezIDs" = entrezIDs)
    rownames(df) <- gpleset@dataTable@table$ID
    return(df)
  }
  
  #look for ensembl ids
  
  if(length(isEnsemBl) != 0){
    
    isEnsemBl <- isEnsemBl[1]
    cat(paste("EnsemBl identifiers found in column", isEnsemBl, ":",
                colnames(gpleset@dataTable@table)[isEnsemBl], "(",
          gpleset@dataTable@columns$Description[isEnsemBl] ,")"), sep = "\n"
        , file = log, append = T)
    
    ensemblIDs <- as.vector(gpleset@dataTable@table[[isEnsemBl]])
    ensemblIDs <- vapply(ensemblIDs, FUN = function(x) ifelse(nchar(x) == 11, paste("ENSG", x, sep = ""), x), FUN.VALUE = character(1), USE.NAMES = FALSE)
    entrezIDs <- getEntrezIDS(genebank, as.list(org.Hs.egACCNUM2EG))
    df <- data.frame("entrezIDs" = entrezIDs)
    rownames(df) <- gpleset@dataTable@table$ID
    return(df)
  }
  isGeneSymbol <- grep("gene symbol" , gpleset@dataTable@columns$Description, ignore.case=TRUE)
  if(length(isGeneSymbol) != 0){
  
    isGeneSymbol <- isGeneSymbol[1]

     cat(paste("Gene symbol identifiers identifiers found in column", isGeneSymbol, ":",
                colnames(gpleset@dataTable@table)[isGeneSymbol], "(",
          gpleset@dataTable@columns$Description[isGeneSymbol] ,")"), sep = "\n"
        , file = log, append = T)

    geneSymbol <- as.vector(gpleset@dataTable@table[[isGeneSymbol]])
    entrezIDs <- getEntrezIDS(geneSymbol, as.list(org.Hs.egALIAS2EG))
    df <- data.frame("entrezIDs" = entrezIDs)
    rownames(df) <- gpleset@dataTable@table$ID
    return(df)
  
  }

  cat("No annotation columns found for genebank, geneid or ensembl", sep = "\n", file = log, append = T)
}


idConversion <- function(element){
  
  # for probes that have more than 1 id (element are the ids)
  # return NA. Its a hot topic what to do with them.
  # For this automatic pipeline we assume we simply cannot know
  # for certain which gene is being expressed and remove them
  #transform null to na, na are more manageble
  
  if(is.null(element)){
    return(NA)
  }
  
  if (length(element) == 1 & any(!is.na(element))){
    #any, while not neccessarily strict due to AND operator, is more stable
    return(element)
  }
  return(NA)
}


getEntrezIDS <- function(identifiers, mapping){
  

  EntrezIds <- as.list(mapping)
  toEntrezIds <- EntrezIds[identifiers]
  entrezIdsOne <- unlist(lapply(toEntrezIds, function(x) idConversion(x)))
  return(entrezIdsOne)
 
}


getMaxMean <- function(x){
  
  #x is a dataframe of sample and entrezIDs columns 
  #if multiple rows correspond to a single gene, return the one
  #with highest mean
 
  idList <-  lapply(split(x[, colnames(x) != "entrezIDs", drop = F], f = x$entrezIDs), FUN=function(x) return(names(which.max((abs(rowMeans(x)))) )))

  return(x[unlist(idList), ])
}

create_dir <- function(name){
  #checks to see if name is a directory, if not
  #creates it. If it exists and is not a directory
  #tough luck
  if(!dir.exists(name)){
    dir.create(name)
  }
}

getBINGSpaceGenes <- function(signGenes, col = "entrezIDs",probeToEntrezIDs = F,topGenes = 150, name, log = "", direc = "."){
  
  #probeToEntrezIDs: dataframe of rows:probeset IDs, one column: entrez gene IDs
  #signGenes: dataframe of rows: probeset IDs, columns logFC
  #generate gmt fileto use in clue io
  #clue io uses only a subset of genespace, as well as a limit of 150 genes
  #subset those allowed, and if there are more than 150, sort based on logFC
  if(probeToEntrezIDs){
    upgenes <-  merge(subset(signGenes, signGenes$logFC > 0), probeToEntrezIDs , 0)
    rownames(upgenes) <- upgenes$Row.names
    upgenes$Row.names <- NULL
    
    dngenes <-  merge(subset(signGenes, signGenes$logFC < 0), gseset@featureData@data, 0)
    rownames(dngenes) <- dngenes$Row.names
    dngenes$Row.names <- NULL
  }
  else{
    signGenes[, col] <- as.integer(signGenes[, col]) 
    upgenes <- subset(signGenes, signGenes$logFC > 0)
    dngenes <- subset(signGenes, signGenes$logFC < 0)
  }
  
  
  upgenes <- upgenes[!is.na(upgenes[, col, drop = F]), ]
  upgenes <- upgenes[upgenes[, col] %in% BINGspace$Entrez.ID, ]
  
  if (nrow(upgenes) > 150){
    upgenes <- upgenes[order(-1*upgenes$logFC), col][1:150]
    
  }else if(nrow(upgenes) < 10){
    
    cat(paste("Less than 10 significantly differentially expressed upgenes found in BINGspace for ", name, ". Analysis stopped", sep = ""), sep = "\n"
        , file = log, append = T)
    return(NULL)
    
  }else{
    
    upgenes <- upgenes[order(-1*upgenes$logFC), col]
  }
  
  
  dngenes <- dngenes[!is.na(dngenes[, col, drop = F]), ]
  dngenes <- dngenes[dngenes[, col] %in% BINGspace$Entrez.ID, ]
  
  writeDngenes <- TRUE
  if (nrow(dngenes) > 150){
    dngenes <- dngenes[order(dngenes$logFC), col][1:150]
    
    
  }else if(nrow(dngenes) < 10){
    
    cat(paste("Less than 10 significantly differentially expressed dngenes found in BINGspace for ", name, ". Skipping dngenes list", sep = ""), sep = "\n"
        , file = log, append = T)
    writeDngenes <- FALSE
    
  }else{
    
    dngenes <- dngenes[order(dngenes$logFC), col]

  }
  
  upgenes <- c(paste(name, "_up", sep = ""), "TAG", upgenes)
  dngenes <- c(paste(name, "_dn", sep = ""), "TAG", dngenes)
  
  fwrite(as.list(upgenes), file = paste(direc,"/gmtFiles/",name, "UP.gmt", sep = ""), sep = "\t")
  
  if(writeDngenes){
    fwrite(as.list(dngenes), file = paste(direc,"/gmtFiles/",name, "DN.gmt", sep = ""), sep = "\t")
    cat(paste("_up and _dn GMT files for ",direc,"/gmtFiles/", name, " created.", sep = ""), sep = "\n"
        , file = log, append = T)
  }
  else{
    
    cat(paste("_up GMT files for ", name, " created.", sep = ""), sep = "\n"
        , file = log, append = T)
  }
}



generalizedTtest<- function(annotatedGseset, designMatrix, contMatrix, id_matrix, log = "", direc, file ){
  
  designMatrix <- designMatrix[order(match(paste("GSM", rownames(designMatrix), sep = ""), colnames(annotatedGseset))), , drop = F]
  #reorder the design matrix so row ordering matches column ordering of annotatedGseset
  #limma will not look for rownames
  exprs(annotatedGseset) <- log2(exprs(annotatedGseset))
  fit <- lmFit(annotatedGseset, designMatrix)
  fit<-contrasts.fit(fit, contMatrix)
  fit<-eBayes(fit)
  
  for(coeff in colnames(contMatrix)){
    
    GSMs <- names(which(designMatrix[,coeff] != 0))
    node <- id_matrix[ id_matrix$GSM %in% strtoi(GSMs), "node_id"  ][1]
    tt <-topTable(fit, coef = coeff, number= nrow(annotatedGseset), adjust.method="BH", sort.by="logFC")
    summarizedGenes <- getMaxMean(tt[ , c("entrezIDs", "logFC"), drop = F])
    
    write.table(tt, file = paste( direc,"/topTable/",  node, "___", file, "_top_table_", coeff ,sep = "") )
    getBINGSpaceGenes(signGenes = summarizedGenes, probeToEntrezIDs = F, name = paste(node, "_",file,coeff, sep = "_"),log = log, direc = direc )
    
  }
}

