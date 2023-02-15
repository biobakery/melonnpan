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

############################################################################################
# Function borrowed from: https://github.com/MRCIEU/metaboprep/blob/master/R/rntransform.R #
############################################################################################

#' Rank-based inverse normal tranformation
#'
#' This function rank normal transforms a vector of data. The procedure is built off of that provided in the GenABEL pacakge.
#'
#' @param y a numeric vector which will be rank normal transformed
#' @param split_ties a binary string of FALSE (default) or TRUE indicating if tied values, of the same rank, should be randomly split giving them unique ranks.
#'
#' @keywords rank normal transformation
#'
#' @importFrom stats qnorm
#'
#' @return returns a numeric vector, with the same length as y, of rank normal transformed values
#'
#' @export
#'
#' @examples
#' ## simulate a negative binomial distribution of values
#' nb_data = rnbinom(500, mu = 40, size = 100)
#' ## rank normal transform those values
#' rnt_data = rntransform( nb_data , split_ties = TRUE )
#'
rntransform = function( y, split_ties = FALSE ){
  ## verify that x is a vector.
  if(is.vector(y) == FALSE){ stop("id_outliers() parameter y is expected to be a vector of data.") }
  
  ## verify some variablity
  if (length(unique(y)) == 1){ stop("trait is monomorphic") }
  
  if (length(unique(y)) == 2){ stop("trait is binary") }
  
  ## identify and remove the NAs
  w = (!is.na(y))
  x <- y[w]
  
  ## empty vector of Z values
  z = rep(NA, length(w))
  
  ## z-transform
  temp <- (x - mean(x))/sd(x)
  
  ## define z values for all observations including the NAs
  z[w] = temp
  
  ## rank the data
  if(split_ties == TRUE){
    rnt <- rank(z, ties.method = "random") - 0.5
  } else {
    rnt <- rank(z) - 0.5
  }
  ## insure the NAs remain so
  rnt[is.na(z)] <- NA
  
  ## inverse
  rnt <- rnt/(max(rnt, na.rm = T) + 0.5)
  ## quantile normalize
  rnt <- stats::qnorm(rnt)
  
  return(rnt)
}
