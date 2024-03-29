#' Model-based Genomically Informed High-dimensional Predictor 
#' of Microbial Community Metabolite Profiles
#'
#' Predict metabolites from microbial sequence features' abundances. 
#' Both measurements are expected to be normalized before using MelonnPan.
#' @param metab Training data of metabolite relative abundances (normalized). 
#' Must have the exact same rows (subjects/samples) as metag. 
#' @param metag Training data of metagenomics sequence features' relative abundances (normalized).
#' Must have the exact same rows (subjects/samples) as metab. 
#' @param output The output folder to write results.
#' @param alpha Grid of alpha values between 0 and 1. Default is `seq(0.05, 0.95, 0.05)`. 
#' @param lambda.choice Choice of optimal lambda ('lambda.min' or 'lambda.1se'). Default is 'lambda.1se'.  
#' @param nfolds Number of folds for internal cross-validation. Default is 10.
#' @param correction Multiplicity adjustment method, same as 'p.adjust'. Default is 'fdr'.
#' @param method Method to correlate measured and predicted metabolites ('spearman' or 'pearson'). Default is 'spearman'.
#' @param cores Number of cores to use for parallel processing. Default is 4.
#' @param seed Specify the arbitrary seed value for reproducibility. Default is 1234.
#' @param cutoff Q-value threshold for significant prediction. Default is 0.05.
#' @param verbose Should detalied message be printed. Default is TRUE.
#' @param no.transform.metab Should arcsine square root transformation (AST) be turned off for `metab`? 
#' Default is FALSE. If FALSE, `metab` must be proportional data ranging from 0 to 1.
#' @param no.transform.metag Should rank-based inverse normal (RIN) transformation be turned off for `metag`? 
#' Default is FALSE. 
#' @param discard.poor.predictions Should predictions with model size = 1 be discarded? Default is FALSE.
#' @param plot Should CV error as a function of lambda be plotted. Default is FALSE.
#' @param outputString Names of the three output files. Default is 
#' 'MelonnPan_Training_Summary', 'MelonnPan_Trained_Weights', and 'MelonnPan_Trained_Metabolites'.
#' @keywords metabolite prediction, microbiome, metagenomics, elastic net, metabolomics
#' @export
melonnpan.train<-function(
              metab, 
              metag, 
              output,
              alpha = seq(0.05, 0.95, 0.05),
              lambda.choice = 'lambda.1se',
              nfolds = 10,
              correction = 'fdr',
              method = 'spearman',
              cores = 4,
              seed = 1234,
              cutoff = 0.05,
              verbose = TRUE,
              no.transform.metab = FALSE,
              no.transform.metag = FALSE,
              discard.poor.predictions = FALSE,
              plot = FALSE,
              outputString = c("MelonnPan_Training_Summary", 
                             "MelonnPan_Trained_Weights", 
                             "MelonnPan_Trained_Metabolites")){
  
  ##############################
  # Read in the input data (X) #
  ##############################
  # if a character string then this is a file name, else it 
  # is a data frame
  if (is.character(metag)){
    temp.metag <-data.frame(readTable(metag))
  } else{    
    temp.metag<-metag
  }
  
  ##############################
  # Read in the input data (Y) #
  ##############################
  # if a character string then this is a file name, else it 
  # is a data frame
  if (is.character(metab)){
    temp.metab <-data.frame(readTable(metab))
  } else{    
    temp.metab<-metab
  }
  
  
  # Re-assign names
  metab<-temp.metab
  metag<-temp.metag
  
  # Dimensionality Check
    if(all(rownames(metab)==rownames(metag))==FALSE)
      stop("Both training datasets should have the same rownames.")

  # Proportionality Check
    # if(any(metab<0)||any(metab>1)|| any(metag<0)||any(metag>1))
    #   stop("All measurements should be normalized to proportions.")
  
  # Create an output folder and figures folder if it does not exist
  if (!file.exists(output)) {
    print("Creating output folder")
    dir.create(output)
  }
   
  # Setting Up Parallel Processing
  cl <- makeCluster(cores)
  registerDoParallel(cl, cores = cores)
  set.seed(seed) # This is for reproducibility
  
  # Use Pre-computed Folds # THIS ENSURES THAT THE RESULTS ARE NOT RANDOM
  foldid <- sample(rep(seq(nfolds), length.out = nrow(metab)))
  
  # Sequence features' Transformation (DEFAULT Rank-Based Inverse Normal)
  if (!no.transform.metag) metag<-apply(metag, 2, melonnpan::rntransform)

  ###############################################################
  # Run Elastic Net (Default Lambda and Custom Alpha Sequence ) #
  ###############################################################
  
  # Metabolite Transformation
  if (!no.transform.metab) {
    if(any(metab<0)||any(metab>1)){
      stop("All metabolite measurements should be normalized to proportions when AST transformation is desired")
    }
    metab<-apply(metab, 2, melonnpan:::ArcSin)
  }
  
  # Elastic Net Regularization
  DD<-CV.ENET(metab=metab, metag=metag, alpha = alpha, lambda.choice = lambda.choice, nfolds=nfolds, foldid = foldid, verbose=verbose, plot=plot, outputDirectory=output)
  
  # Back Transformation to Project Predicted Values in the Original Space
  if (!no.transform.metab){
    pred<-apply(DD$pred, 2, melonnpan:::SqSin)
  } else{
    pred<-DD$pred
    }
  
  # Summarize Results
  weight.matrix<-DD$weight.matrix
  colnames(weight.matrix)<-colnames(metab)
  rownames(weight.matrix)<-c('Intercept', colnames(metag))
  compounds<-melonnpan.summarize(metab=metab, pred=pred, method=method, correction=correction, cutoff=cutoff)
   
  # Remove IDs from Weight and Predicted Matrices If Not Present in Compounds
  weight.matrix<-weight.matrix[,intersect(colnames(weight.matrix),  compounds$ID)]
  pred<-pred[,intersect(colnames(pred), compounds$ID)]
   
  # Save model size
  modelSize<-as.data.frame(colSums(weight.matrix!=0)) - 1
  modelSize_DF<-rownames_to_column(modelSize, 'ID')
  colnames(modelSize_DF)<-c('ID', 'Model.Size')
  
  # Add structure to weight and predicted matrix
  final_weight_matrix<-as.data.frame(weight.matrix)
  final_weight_matrix<-rownames_to_column(final_weight_matrix, 'ID')
  final_predicted_matrix<-as.data.frame(pred)
  final_predicted_matrix<-rownames_to_column(final_predicted_matrix, 'ID')
  
  if (is.null(compounds)){
    cat("MelonnPan was unable to find well-predicted metabolites. No output returned.")
    } 
  else{
    
    # Combine with Model Size
    final_compounds_table<-merge(compounds, modelSize_DF, 'ID')
    
    if(discard.poor.predictions){
      
      ############################
      # Subset valid predictions #
      ############################
      
      # Remove metabolites with model size <=1
      final_compounds_table<-final_compounds_table[final_compounds_table$Model.Size>1,]
      final_compounds_table$Q.value<-p.adjust(final_compounds_table$P.value)
      final_compounds_table<-final_compounds_table[order(final_compounds_table$Q.value, decreasing = FALSE),]
      
      ################################################
      # Subset weight matrix and predicted compounds #
      ################################################
      
      final_predicted_matrix<-column_to_rownames(final_predicted_matrix, 'ID')
      final_predicted_matrix<-final_predicted_matrix[, final_compounds_table$ID]
      final_predicted_matrix<-rownames_to_column(final_predicted_matrix, 'ID')
      
      final_weight_matrix<-column_to_rownames(final_weight_matrix, 'ID')
      final_weight_matrix<-final_weight_matrix[, final_compounds_table$ID]
      final_weight_matrix<-rownames_to_column(final_weight_matrix, 'ID')
    }
    
    # Save Outputs
    write.table(final_compounds_table,
                  paste(output, paste(outputString[[1]], '.txt', sep=''), sep='/'),
                  row.names=F,col.names=T,quote=F,sep="\t")
    write.table(final_weight_matrix,
                paste(output, paste(outputString[[2]], '.txt', sep=''), sep='/'),
                row.names=F,col.names=T,quote=F,sep="\t")
    write.table(final_predicted_matrix,
                paste(output, paste(outputString[[3]], '.txt', sep=''), sep='/'),
                row.names=F,col.names=T,quote=F,sep="\t")
    cat("The output tables are in", output)
    }
}


