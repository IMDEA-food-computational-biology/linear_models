Sys.setenv("VROOM_CONNECTION_SIZE" = 131072 * 8)
source("detectDesign.R")
source("normalizedArraysFunctions.R")
source("HTS_functions.R")

args = commandArgs(trailingOnly=TRUE)
DM <- args[1]
data_matrix <- args[2]
nodes <- args[3]

#filename suffix for outputs
f <- strsplit(data_matrix, "/")[[1]]
f <- f[length(f)]
f <- strsplit(f, "_")[[1]][1]

con <- file(DM,"r")
first_line <- readLines(con,n=1)

ret <- detect_design(DM, log = "", sep = "\t")
#print(ret$designMatrix)
print(ret$contMatrix)
#checker wheter HTS or microarrays
if(first_line == "#HTS"){

	counts <- prepare_matrix(data_matrix, sep = "\t", header = 1) 
	annotatedEset <- proccessHTS(counts, ret$designMatrix)
}else{
	annotatedEset <- annotateGSE(data_matrix)
}


node_ids <- read.table(nodes, sep = "\t", header = T)
generalizedTtest(annotatedEset, ret$designMatrix, ret$contMatrix,node_ids,"",".", f)
