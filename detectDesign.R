# each design has its own function
# lots of code repetition for the special cases
# generateSpecialDesign.R has one for all the special types
# metaDataFile = "~/Marco_TFM/Marco_TFM/Data/metaData/GSE111073sample_info.txt"
# isPaired = T
# intercept = F
#install.packages("stringr")
library(stringr)

getMeta <- function(meta,...){

  "
  If meta is string, read table
  If dataframe, just return it
  "
  if(is.character(meta) & length(meta) == 1){
	return(read.table(meta, header = TRUE, stringsAsFactors = F, ...))
  }else if (is.data.frame(meta)){
  	return(meta)
  }else{

        print("meta must be a file path or data frame")
        return(1)
  }

}

generateFormula <- function(dfName = "", varnames, intercept = F){
  
  if(intercept){
    formula <- "~1+"
  }else{
    formula <- "~0+"
  }
  
  i = 0
  
  for(var in varnames){
    
    if( i == 0){
      formula <- paste(formula,  dfName,"$", var,sep = "")
    }else{
      formula <- paste(formula, " + ",dfName,"$", var, sep = "")
    }
    i <- i + 1
    
  }
  return(as.formula(formula))
}

two_color_coeffs <- function(gsm, m){

    m <- m[m[,"varDataG$GSM"] == gsm,]
    return(m[m[,"varDataG$colorCy5"] == 1,,drop = F] - m[m[,"varDataG$colorCy5"] == 0,,drop = F])
}

build_two_color <- function(designMatrix, f){

    gsms <- unique(designMatrix[, "varDataG$GSM"])
    newDM <- do.call(rbind.data.frame, lapply(gsms,FUN = two_color_coeffs, m = designMatrix))

    newDM <- newDM[, -grep("^varDataG\\$(GSM|colorCy5)", colnames(designMatrix))]
    rownames(newDM) <- gsms

    rank <- qr(newDM)$rank
    if(rank != (ncol(newDM))){

        cat("WARNING: rank of design is one less than full rank. Assuming controls are redundant in 2 color experiment\n", file = f, append = T)
        #print(colnames(newDM))
        newDM <- newDM[, -grep( "^varDataG\\$treatmentcontrol" , colnames(newDM)), drop = F]
    }
    return(as.matrix(newDM))
}

get_rows_and_color <- function(designMatrix, two_color, varData, file = ""){


  "
  If one color, assign rownames to GSM, remove GSM columns
  If two color, create DM by subtracting cy3 from cy5
  "
  if(two_color){

      designMatrix <- build_two_color(designMatrix, file)
  }else{

          designMatrix <- designMatrix[, -grep("^varDataG\\$GSM", colnames(designMatrix)), drop = F]
          rownames(designMatrix) <- varData$GSM
  }
  return(designMatrix)
}

generateStandardDesignMatrix <- function(meta, isPaired = F, intercept = F, two_color = F,file = "",...){

  varData <- getMeta(meta,...)
  if(isPaired | "Id" %in% colnames(varData)){
    varData$Id <- as.factor(varData$Id)
  }

  n <- grep("treatment", colnames(varData))
  assign("varDataG", varData, envir = .GlobalEnv) #needed because model.matrix looks in the global environment
  designMatrix <- model.matrix( generateFormula("varDataG", colnames(varData), intercept = intercept))
  designMatrix <- get_rows_and_color(designMatrix, two_color, varData)

  colnames(designMatrix) <- vapply(strsplit(colnames(designMatrix), "\\$" ), function(x) make.names(x[2]), FUN.VALUE = character(1))
  return(designMatrix)
}



generateStandardContrastMatrix <- function(designMatrix,file = "", ...){

  treatmentNames <- colnames(designMatrix)[grep("^treatment.*", colnames(designMatrix), ignore.case = T)]
  t <- treatmentNames[treatmentNames != "treatmentcontrol"]
  contMatrix <- matrix(0, nrow = ncol(designMatrix), ncol = length(t),dimnames = list(colnames(designMatrix), t) )
  control <- which(colnames(designMatrix) %in% c("treatmentcontrol"))
  if(length(control) == 0){

      cat("WARNING: No control found, check if control has been made redundant in 2 color experiment\n", file = file, append = T)
  }
  contMatrix[control, ] <- -1

  for(i in t){

    contMatrix[i, i] <- 1

    }

  return(contMatrix)
}

# diff time points, each w/a control

generateDesignTimePoint <- function(meta, isPaired = F, intercept = F, varName= "time point", two_color = F,file = "",...){


  ## diff time points
  ## no 0 time point
  ## each time point has its own control or general control
  varData <- getMeta(meta,...)
  "  
  varData <- data.frame(vapply(X=varData, FUN = function(x) make.names(x), FUN.VALUE = character(nrow(varData)))
                        ,row.names = row.names(varData))
  varName <- make.names(varName)"
  varData$treatment <- str_replace_all(varData$treatment, pattern = "(\\.| |/)", "_")

  if(isPaired | "Id" %in% colnames(varData)){
    varData$Id <- as.factor(varData$Id)
  }
  treatmentVar <- paste(varData$treatment, varData[ , varName], sep = "___")
  varData <- data.frame( treatment = treatmentVar,  varData[,!colnames(varData)  %in% c("treatment", varName), drop = F],
                            stringsAsFactors = F)


  assign("varDataG", varData, envir = .GlobalEnv)
  designMatrix <- model.matrix( generateFormula("varDataG", colnames(varData), intercept = F) )
  designMatrix <- get_rows_and_color(designMatrix, two_color,  varData, file = file)
  colnames(designMatrix) <- vapply(strsplit(colnames(designMatrix), "\\$" ), function(x) make.names(x[2]), FUN.VALUE = character(1))
  return(designMatrix)
  }

generateContrastTimePoint <- function(designMatrix, varName = "time point",file = "",...){

  treatment_indexes <- grep(x = colnames(designMatrix),pattern =  "^treatment(?!control)", perl = T)
  treatment_varNames <- strsplit(colnames(designMatrix)[treatment_indexes], split = "___" )

  treatment <- unique(lapply(treatment_varNames, function(x) x[[1]]))
  timepoints <- unique(lapply(treatment_varNames, function(x) paste(x[2:length(x)], collapse = "___")))

  contMatrix <- matrix( data = 0, nrow = ncol(designMatrix), ncol = length(treatment)*length(timepoints))
  colnames(contMatrix) <- paste(rep(treatment, each = length(timepoints)), timepoints,sep = "___")
  rownames(contMatrix) <- colnames(designMatrix)
  for(co in colnames(contMatrix)){
    col_split <- strsplit(co, "___")[[1]]
    tp <- paste(col_split[2:length(col_split)], collapse = "___")
    control <- which(colnames(designMatrix) %in% paste("treatmentcontrol", tp, sep = "___"))
    print(paste("treatmentcontrol", tp, sep = "___"))
    print(co)
    print("::")
    contMatrix[co, co] <- 1
    if(length(control) == 0){

        cat("WARNING: No control found, check if control has been made redundant in 2 color experiment\n", file = file, append = T)
    }
    contMatrix[control, co] <- -1

  }

  return(contMatrix)

}



generateDesign0asControl <- function(meta, isPaired = F, varName = "time_point", two_color= F,...){

  #diff time points/conc, 0 is control
  #control can be shared across treatments or not
  varData <- getMeta(meta,...)
  #print(varName)
  #print(make.names(varName))
  "varData <- data.frame(vapply(X=varData, FUN = function(x) make.names(x),
                               FUN.VALUE = character(nrow(varData)))
                        ,row.names = row.names(varData))
  #print(colnames(varData))
  varName <- make.names(varName)"
  varData$treatment <- str_replace_all(varData$treatment, pattern = "(\\.| |/)", "_")

  if(isPaired | "Id" %in% colnames(varData)){
    varData$Id <- as.factor(varData$Id)
  }
  #print(varData)
  #print(colnames(varData))
  #print(varName)
  #print(varData[ , varName])
  treatmentVar <- paste(varData$treatment, varData[ , varName], sep = ".")
  varData <- data.frame( treatment = treatmentVar,  varData[,!colnames(varData)  %in% c("treatment", varName), drop = F],
                           stringsAsFactors = F)

  assign("varDataG", varData, envir = .GlobalEnv)
  designMatrix <- model.matrix( generateFormula("varDataG", colnames(varData), intercept = F) )
  designMatrix <- get_rows_and_color(designMatrix, two_color, varData)
  colnames(designMatrix) <- vapply(strsplit(colnames(designMatrix), "\\$" ), function(x) make.names(x[2]), FUN.VALUE = character(1))
  #rownames(designMatrix) <- rownames(varData)
  #print(unique(treatmentVar))
  return(designMatrix)
  }

generateContrast0asControl  <- function(designMatrix, ...){

  treatment_indexes <- grep(x = colnames(designMatrix),pattern =  "^treatment", perl = T)
  treatment_varNames <- strsplit(colnames(designMatrix)[treatment_indexes], split = "\\." )

  treatmentCompounds <- unique(unlist(lapply(treatment_varNames, function(x) x[1])))
  treatmentVars <- unique(unlist(lapply(treatment_varNames, function(x) paste(x[2:length(x)], collapse = "."))))

  if("treatmentcontrol" %in% treatmentCompounds){

	  treatmentCompounds <- treatmentCompounds[treatmentCompounds != "treatmentcontrol"]
	  treatmentVars <- treatmentVars[treatmentVars != "X0"]
	  control_idxs <- grep(x = colnames(designMatrix),pattern =  "^treatmentcontrol", perl = T)
	  contrast_idxs <- setdiff(treatment_indexes, control_idxs)
	  contrast_columns <- colnames(designMatrix)[contrast_idxs]
	  contMatrix <- matrix( data = 0, nrow = ncol(designMatrix),
                           ncol = length(contrast_columns))

	  colnames(contMatrix) <- contrast_columns
          rownames(contMatrix) <- colnames(designMatrix) #already has
	  contMatrix["treatmentcontrol.0", ] <- -1
	  for(co in colnames(contMatrix)){
   	 	contMatrix[co, co] <- 1
	 }
  }else{

  	#diff conc for each compound: redo for everything!
    	treatmentVar <- unique(colnames(designMatrix)[treatment_indexes])
    	contMatrix <- matrix( data = 0, nrow = ncol(designMatrix), ncol = length(treatmentVar) - length(treatmentCompounds))
  	colnames(contMatrix) <- treatmentVar[-grep(pattern = "\\.0$", x = treatmentVar)]
 	rownames(contMatrix) <- colnames(designMatrix) #already has
  	for(co in colnames(contMatrix)){

		control <- sub("\\..*", ".0", co)
		contMatrix[co, co] <- 1
		contMatrix[control, co] <- -1
  	}
  }

  return(contMatrix)
}

generateDesignBeforeAfter <- function(meta, isPaired = F, varName = "time point", two_color = F,...){

  #before after studies
  varData <- getMeta(meta,...)
  "varData <- data.frame(vapply(X=varData, FUN = function(x) make.names(x), FUN.VALUE = character(nrow(varData)))
                        ,row.names = row.names(varData), stringsAsFactors = F)

  #in case some treatment has a dot due to bad naming, change to _, otherwise will interfere
  #with value assingmen in contMatrix
  varName <- make.names(varName)"
  varData$treatment <- str_replace_all(varData$treatment, pattern = "(\\. | |/|-)", "_")


  if(isPaired | "Id" %in% colnames(varData)){
    varData$Id <- as.factor(varData$Id)
  }

  treatmentVar <- paste(varData$treatment, varData[ , varName], sep = "___")
  varData <- data.frame( treatment = treatmentVar,  varData[,!colnames(varData)  %in% c("treatment", varName), drop = F],
                            stringsAsFactors = F)

  assign("varDataG", varData, envir = .GlobalEnv)
  designMatrix <- model.matrix( generateFormula("varDataG", colnames(varData), intercept = F) )
  designMatrix <- get_rows_and_color(designMatrix, two_color, varData)
  colnames(designMatrix) <- vapply(strsplit(colnames(designMatrix), "\\$" ), function(x) make.names(x[2]), FUN.VALUE = character(1))
  #rownames(designMatrix) <- rownames(varData)
  return(designMatrix)
  }

generateContrastBeforeAfter <- function(designMatrix, ...){


  treatment_indexes <- grep(x = colnames(designMatrix),pattern =  "^treatment", perl = T)
  treatment_varNames <- colnames(designMatrix)[treatment_indexes]
  treatment_varNames_L <- strsplit(treatment_varNames, split = "___" )

  #treatment <- unique(lapply(treatment_varNames, function(x) x[[1]]))
  timepoints <- unique(lapply(treatment_varNames_L, function(x) paste(x[2:length(x)], collapse = "___")))

  treatment <- unique(treatment_varNames)
  treatment <- treatment[-c(grep(pattern = "^treatmentcontrol", x = treatment),
                            grep(pattern = "___0$", x = treatment))]

  timepoints <- timepoints[which(timepoints != "X0")]

  contMatrix <- matrix( data = 0, nrow = ncol(designMatrix), ncol = length(treatment))
  colnames(contMatrix) <- treatment
  rownames(contMatrix) <- colnames(designMatrix)

  contMatrix["treatmentcontrol___0", ] <- 1 #here control means placebo
  for(co in colnames(contMatrix)){

    col_split <- strsplit(co, "___")[[1]]
    tr <- col_split[1]
    tp <- paste(col_split[2:length(col_split)], collapse = "___")

    control <- paste(tr, "___0", sep = "") #here control refers to time point 0
    cat("co", co, "\n")
    contMatrix[co, co] <- 1

    if(control %in% row.names(contMatrix)){
	#sometimes we dont have a time point 0 for a given condition
	#which os not ideal
	#but here we just treat those case as if it were a diff time point
	#study, with the benefit that control has its time
	# effect 0 subtracted 
	#TODO: what to do with control 0:subtarct? nothing? for now nothing
	#(remove control 0 coefficient)
	cat("control", control, "\n")
    	contMatrix[control, co] <- -1
    }else{
    	contMatrix["treatmentcontrol___0", co] <- 0
    }
    cat("tp", tp, "\n")
    contMatrix[paste("treatmentcontrol", tp, sep = "___"), co] <- -1
  }


  return(contMatrix)
}

detect_design <- function(metaDataFile, log = "",...){

	cat("^", metaDataFile, "\n", file = log, append = T)
	con <- file(metaDataFile,"r")
	first_line <- readLines(con,n=1)
	close(con)
	varData <- read.table(metaDataFile, header = TRUE, ...)
	cols <- colnames(varData)

	admitted_lines <- c("#HTS","#One color array","#Two channel array")
	if(!(first_line %in% admitted_lines)){

		cat("First line must be one of these three",admitted_lines,"\n",file = log, append = T)
		return(1)
	}

	tp <- grep(pattern = "^time", cols)
	conc <- grep(pattern = "^concentration", cols)
	
	cat("Dealing with ", first_line, "\n", file = log, append = T)
	two_channel = first_line == "#Two channel array"
	
	if(length(tp) == 0 & length(conc) == 0){

	    cat("Simple study\n", file = log, append = T)
	    DM <- generateStandardDesignMatrix(varData,two_color = two_channel,...)
	    cont <- generateStandardContrastMatrix(DM)
	}else{
	    if(two_channel){
		'
		cat("No implementation yet for non standard two channel studies\n", file = log, append = T)
            	return(1)'
	    }
	    if( length(tp) > 0 & length(conc) > 0 ){
		    cat("Concentration and time detected. Concatenating treatment and concentration. If this is undesired, your study might not be of accepted type.\n", file = log, append = T)
		    varData$treatment <- paste(varData$treatment, varData$concentration, sep = "_")
		    varData[grep(pattern = "^control", x = varData$treatment), "treatment"] <- "control"
		    varData <- varData[, !(colnames(varData) %in% c("concentration"))]
		    conc <- integer(0)
		}
	    column <- ifelse(length(tp), tp,conc)
	    varName <- cols[column]
	    if (!("0" %in% varData[, column])){

		cat("Each time point its own control\n", file = log, append = T)
		DM <- generateDesignTimePoint(varData, varName = varName, two_color = two_channel,...)
		cont <- generateContrastTimePoint(DM)
	    }else if (!("control" %in% varData$treatment)){

		cat("One control per time and treatment\n", file = log, append = T)
		DM <- generateDesign0asControl(varData,varName = varName, two_color = two_channel,...)
		cont <- generateContrast0asControl(DM)
	    }else{

		controls <- which(varData$treatment == "control")
		zeros <- which(varData[, column] == "0")

		 if( (length(controls) == length(zeros)) && all(controls == zeros)){

		    cat("One control per treatment for all time points\n", file = log, append = T)
		    DM <- generateDesign0asControl(varData,varName = varName,two_color = two_channel, ...)
		    cont <- generateContrast0asControl(DM)
		}else{

		cat("Before after studies\n", file = log, append = T)
		DM <- generateDesignBeforeAfter(varData,varName = varName,two_color = two_channel,...)
		cont <- generateContrastBeforeAfter(DM)
		}
	}

	}
ret <- list()
ret$designMatrix <- DM
ret$contMatrix <- cont
return(ret)

}
