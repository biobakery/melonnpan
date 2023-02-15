#' Model-based Genomically Informed High-dimensional Predictor 
#' of Microbial Community Metabolite Profiles
#'
#' Predict metabolites from new microbiome samples. 
#' 
#' @param metag Microbial sequence features' relative abundances (matrix) 
#' for which prediction is desired. The sequence features' abundances are expected to 
#' be normalized.
#' @param output The output folder to write results.
#' @param no.transform.metab Should back transformation be turned off for predicted metabolites? 
#' Default is FALSE. If FALSE, it is expected that the proportional data ranging from 0 to 1 was used for training. 
#' @param no.transform.metag Should rank-based inverse normal (RIN) transformation be turned off for `metag`? 
#' Default is FALSE.
#' @param weight.matrix The weight matrix to be used for prediction (optional). 
#' If not provided, by default, a pre-trained weight matrix based on UniRef90 gene 
#' families from the original MelonnPan paper (Mallick et al, 2019) will be used.
#' @param train.metag Quality-controlled training metagenomes against which 
#' similarity is desired (optional). The sequence features' abundances are expected 
#' to be normalized.
#' If not provided, a pre-processed UniRef90 gene family training table from 
#' the original MelonnPan paper (Mallick et al. 2019) will be used.
#' @param criticalpoint A numeric value corresponding to the significance level 
#' to find the top PCs. If the significance level is 0.05, 0.01, 0.005, or 0.001, 
#' the criticalpoint should be set to be 0.9793, 2.0234, 2.4224, or 3.2724, 
#' accordingly. The default is 0.9793 (i.e. 0.05 significance level).
#' @param corr.method Method to correlate new metagenomes and training PCs. 
#' Default is 'pearson'.
#' @keywords metabolite prediction, microbiome, metagenomics, elastic net, metabolomics
#' @export
melonnpan.predict<-function(
              metag,
              output,
              no.transform.metab = FALSE,
              no.transform.metag = FALSE,
              weight.matrix = NULL,
              train.metag = NULL,
              criticalpoint = 0.9793,
              corr.method = 'pearson'){
  
  #################################################
  # Read in the input data (i.e. new metagenomes) #
  #################################################
  # if a character string then this is a file name, else it 
  # is a data frame
  if (is.character(metag)){
    test.metag <-data.frame(readTable(metag))
  } else{    
    test.metag<-data.frame(metag)
  }

  # Sanity check for proportionality 
  # if(any(test.metag<0)||any(test.metag>1))
  #   stop("All measurements should be normalized to proportions.")
 
  # Create an output folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
  
  ##################
  # Calculate RTSI #
  ##################
  
  # Load training metagenomes
  if (is.null(train.metag)){
    train.metag<-melonnpan::melonnpan.training.data
  } else{ 
    # if a character string then this is a file name, else it 
    # is a data frame
    if (is.character(train.metag)){
      train.metag<-data.frame(readTable(train.metag))
    } else{    
      train.metag <-data.frame(train.metag)
    }
  }
  
  # Subset to common IDs
  commonID<-intersect(colnames(train.metag), colnames(test.metag))
  
  # Throw error if no common IDs between training and test data
  if(length(commonID)<1) stop('No common IDs found between training and test data. Execution halted!')
  
  # Common features across datasets
  train<-as.data.frame(train.metag[, commonID])
  test<-as.data.frame(test.metag[, commonID])
  
  # Remove binary features
  ID<-which(colSums(test!=0)>1)
  train<-train[, ID]
  test<-test[,ID]
  rowTrain<-rownames(train)
  rowTest<-rownames(test)
  
  # RIN transformation (For RTSI calculationm original data in still untransformed)
  if (!no.transform.metag) {
    train<-apply(train, 2, melonnpan::rntransform)
    test<-apply(test, 2, melonnpan::rntransform) 
    }
  rownames(train)<-rowTrain
  rownames(test)<-rowTest
  
  # Extract PCs 
  PCA <- prcomp(train)
  
  # Select top PCs based on TW statistic 
  DD<-AssocTests::tw(eigenvalues = PCA$sdev, eigenL = length(PCA$sdev), criticalpoint = criticalpoint) 
  Loadings <- as.data.frame(PCA$rotation[,1:DD$SigntEigenL])
  
  # Calculate pairwise correlation between PCs and samples 
  RTSI<-matrix(nrow=nrow(test), ncol=ncol(Loadings))
  for (i in 1:nrow(test)){
    y<-as.numeric(test[i,])
    for (j in 1:ncol(Loadings)){
      r<-as.numeric(Loadings[,j])
      RTSI[i, j] <- cor(r, y, method = corr.method)
    }
  }
  
  # Add structure to the output
  rownames(RTSI)<-rowTest
  RTSI_Score<-as.data.frame(apply(RTSI, 1, max))
  colnames(RTSI_Score)<-'RTSI'
  RTSI_Score<-tibble::rownames_to_column(RTSI_Score, 'ID')
  write.table(RTSI_Score, 
              file = paste(output, '/MelonnPan_RTSI.txt', sep = ''),
              row.names=F, col.names=T, quote=F, sep="\t") 
  
  #######################
  # Predict metabolites #
  #######################
  
  # Remove binary features (so that RIN transformation is valid)
  retainIDs<-which(colSums(test.metag!=0)>1)
  new.metag<-test.metag[, retainIDs]
  
  # Load weight matrix
  if (is.null(weight.matrix)){
    train.weight<-as.data.frame(t(melonnpan::melonnpan.trained.model))
  } else{ 
    # if a character string then this is a file name, else it 
    # is a data frame
    if (is.character(weight.matrix)){
      train.weight<-as.data.frame(t(readTable(weight.matrix))) # Genes + intercept in columns, metabolites in rows
    } else{    
      train.weight <-as.data.frame(weight.matrix) # Genes + intercept in columns, metabolites in rows
    }
  }

  # Subset by overlapping sequence features
  # i.e. common in both weight matrix and input data (sequence features)
  X<-new.metag[, intersect(colnames(train.weight), colnames(new.metag))]
  X<-X[,which((nrow(X)-colSums(X==0))>1)]
  intercept<-train.weight$Intercept
  test.weight<-train.weight[, intersect(colnames(train.weight), colnames(X))]
  new.weight<-cbind.data.frame('Intercept' = intercept, test.weight)
  compound_names<-rownames(new.weight)
  
  # Apply Rank-based Inverse Normal (RIN) transformation (Separate from RTSI calculation)
  if (!no.transform.metag){
    transf.X<-apply(X, 2, melonnpan::rntransform)
  } else{
    transf.X<-X
  }
  new.X<-cbind('Intercept'=rep(1, nrow(transf.X)), transf.X)
  
  # Carry out prediction in new samples and back-transform if needed
  new.X<-as.matrix(new.X)
  new.weight<-apply(as.matrix.noquote(new.weight),2,as.numeric)
  pred<-new.X%*%t(new.weight)
  
  # Back transformation
  if (!no.transform.metab) pred<-apply(pred, 2, SqSin)
  
  # Add structure to the output
  rownames(pred)<-rownames(new.metag)
  colnames(pred)<-compound_names
  pred<-as.data.frame(pred)
  
  # Remove predictions with constant values
  pred<-pred[, sapply(pred, function(v) var(v, na.rm=TRUE)!=0)] 
  pred<-tibble::rownames_to_column(pred, 'ID')
  write.table(pred, 
              file = paste(output, '/MelonnPan_Predicted_Metabolites.txt', sep = ''),
              row.names=F, col.names=T, quote=F, sep="\t")
  
  # Return 
  return(list(pred = pred, RTSI = RTSI_Score))
}


