#' Plot CV Error As A Function of Lambda
#'
#' This function plots the CV error as a function of lambda from Elastic Net Training.
#' @param CV A glmnet object.
#' @param x Index of glmnet fit.
#' @param alpha Optimal alpha selected.
#' @param outputDirectory Name of the desired output directory. Default is the current directory.
#' @keywords metabolite prediction, microbiome, metagenomics, elastic net, metabolomics
#' 
#' @export
#' 
generateLambdaPlot<-function(CV, x, alpha, outputDirectory){
  pdf(paste(outputDirectory, paste('MelonnPan_LambdaPlot_', x, '.pdf', sep=''), sep=''))
  plot(CV)
  title(paste(x, ', optimal alpha = ', alpha, sep = ''), line = 2.5)
  dev.off()
}

#' Arc Sine Square Root Transformation for Proportional Data
#'
#' This function applies arcsine square root transformation to proportional data.
#' @param x A numerical vector of proportions (must be between 0 and 1).
#' @keywords arcsine, compositional data, proportions, transformation
#' @seealso \code{\link{SqSin}} 
#' @examples
#' ArcSin(runif(100,0,1))
#' 
#' @export
#' 
ArcSin<-function(x){
  if (!is.vector(x) || !(mode(x)=="numeric")){
    stop("'x' must be a numeric vector")
  }
  if (any(x < 0)){
    stop("Some values of x are < 0")
  }
  if (any(x > 1)){
    stop("Some values of x are > 1")
  }
  else{
    return(asin(sqrt(abs(x))))
  }
}



#' Back Transformation of Predicted Values from A Statistical Model with
#' ArcSine Square Root Transformed Compositional Response Variable
#'
#' This function applies back transformation to the predicted values of a statistical model with
#' arcsine square root transformed compositional response variable.
#' @param x A numerical vector of arcsine square root transformed compositions or proportions.
#' @param adjust Should values less than 0 be mapped to 0 and 
#' values greater than pi/2 be mapped to 1? Defaults to TRUE.
#' @keywords arcsine, compositional data, proportions, back transformation
#' @seealso \code{\link{ArcSin}} 
#' @examples
#' 
#' # This is a simple demonstration
#' y <- ArcSin(runif(100,0,1))
#' x <- rnorm(100,0,1)
#' z <- predict(lm(y ~ x))
#' SqSin(z)
#' @export
#' 
SqSin <- function(x, adjust=TRUE){
  if (!is.vector(x) || !(mode(x)=="numeric")){
    stop("'x' must be a numeric vector")
  }
  if(adjust==TRUE){
    return(ifelse(x<0, 0, ifelse(x>pi/2, 1, sin(x)^2)))
  }
  if(adjust==FALSE){
    return(sin(x)^2)
  }
}


#' Read data using a wrapper of data.table
#'
#' This function reads the input data using the data.table package 
#' but returns a matrix with row names (unlike fread). 
#' @param Input Full path of the input data.
#' @keywords metabolite prediction, microbiome, metagenomics, elastic net, metabolomics
#' 
#' @export
#' 
readTable<-function(Input){
  
  # Read Data
  data <- data.table::fread(Input, header=FALSE)
  
  # Remove Dupilcate Row Names 
  rowNames<-data$V1
  data<-dplyr::select(data, -V1)
  
  # Assign Column Names
  colNames<-as.character(unlist(data[1,]))
  colnames(data)<-colNames
  
  # Remove Dupilcate Column Names 
  data<-data[-1,]
  
  # Convert to Matrix
  dat<-as.matrix(as.data.frame(lapply(data, as.numeric)))
  
  # Assign Row Names 
  rownames(dat)<-rowNames[-1]
  
  # Return As Matrix
  return(dat)
}







