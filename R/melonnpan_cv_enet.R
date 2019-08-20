#' Fit a regularized linear model
#'
#' @param metab Training data of metabolite relative abundances. 
#' Should have the exact same rows (subjects/samples) as metag. 
#' @param metag Training data of microbial sequence features' relative abundances.
#' Should have the exact same rows (subjects/samples) as metab. 
#' @param alpha Grid of alpha values between 0 and 1. Default is `seq(0.05, 0.95, 0.05)`. 
#' @param lambda.choice Choice of optimal lambda ('lambda.min' or 'lambda.1se'). Default is 'lambda.1se'.  
#' @param nfolds Number of folds for internal cross-validation. Default is 10.
#' @param foldid A vector of values between 1 and nfold identifying what fold each observation is in. 
#' @param verbose Should progress bar be printed. Default is TRUE.
#' @param plot Should CV error as a function of lambda be plotted. Default is FALSE.
#' @param outputDirectory Name of the desired output directory.
#' @keywords metabolite prediction, mi crobiome, metagenomics, elastic net, metabolomics
#' @export
CV.ENET<-function(metab = metab, 
                metag = metag, 
                alpha = alpha,
                lambda.choice = lambda.choice,
                nfolds = nfolds,
                foldid = foldid,
                verbose=verbose, 
                plot = plot, 
                outputDirectory=outputDirectory){
  
  # Print Progress
  if (verbose==TRUE) {
    cat(c("Job started at:",date()),fill=TRUE)
    cat("\nTraining a MelonnPan model...")
    cat("\nPredicting metabolite compositions using MelonnPan...")
    cat("\nFitting an elastic net regularized linear model...")
    cat("\nFinding the best combination of alpha and lambda...")
    cat('\nk-fold CV:', nfolds)
  }
  
  # Track Time
  start.time <- Sys.time()
  
  # Fit Elastic Net with User-defined Parameters
  tune.cv.enet <- rbind(foreach(r =1:ncol(metab),
                                  .packages=c("glmnet", "melonnpan"),
                                  .combine=rbind,
                                  .inorder =F) %dopar% {
                          cv_run <-cvar.glmnet(x = as.matrix(metag), y = metab[,r], alpha = alpha, nfolds=nfolds, foldid = foldid, lambda.choice = lambda.choice)
                          if (plot==TRUE) generateLambdaPlot(cv_run$opt_cv, colnames(metab)[r], cv_run$opt_alpha, outputDirectory)
                          cv_run_extract <-c('opt_lambda' = cv_run$opt_lambda, 'opt_alpha' = cv_run$opt_alpha)})
  
  # Generate Predictions
  tune.cv.enet<-as.data.frame(tune.cv.enet)
  enet.sol<- sapply(1:ncol(metab), function(r) {glmnet::glmnet(as.matrix(metag),metab[,r], lambda = tune.cv.enet$opt_lambda[r], alpha = tune.cv.enet$opt_alpha[r])})
  enet.coeff<-sapply(1:ncol(metab), function(r){glmnet::predict.glmnet(enet.sol[,r], type="coefficients", s = tune.cv.enet$opt_lambda[r])})
  weight.matrix<-data.matrix(Reduce(cbind, enet.coeff))
  pred<-cbind(rep(1, nrow(metag)), as.matrix(metag))%*%as.matrix(weight.matrix)
  colnames(pred)<-colnames(metab)
  rownames(pred)<-rownames(metab)
  stop.time <- Sys.time()
  time<-round(difftime(stop.time, start.time, units="min"),3)
  if (verbose==TRUE) cat(c("\nJob finished at:",date()),fill=TRUE)
  if (verbose==TRUE) cat("Computational time:", time, "minutes \n")
  return(list(pred=pred, weight.matrix=weight.matrix))
}


###############################################################
# Elastic Net CV with Simultaneous Tuning of Alpha and Lambda #
###############################################################

cvar.glmnet<-function(x, y, nfolds, foldid, alpha, lambda.choice)
{
  a1<-lapply(alpha, .cvfunc, xmat=x, ymat =y, nfolds = nfolds, foldid = foldid, parallel = TRUE)
  DD <- plyr::ldply(a1, extractGlmnetInfo)
  DD$alpha<-alpha
  if (lambda.choice == 'lambda.1se') {
    min_index<-which(DD$cv.1se==min(DD$cv.1se))
    }
  else {
    min_index<-which(DD$cv.1se==min(DD$cv.1se))
    }
  opt_alpha<-DD$alpha[min_index]
  opt_cv<-a1[[min_index]]
  opt_lambda<-opt_cv$lambda.1se
  return(list(opt_cv = opt_cv, opt_lambda = opt_lambda, opt_alpha = opt_alpha))
}  

##################################################
# Perform Elastic Net CV for Each Value of Alpha #
##################################################

.cvfunc <- function(a, xmat, ymat, nfolds, foldid, parallel, ...)
{
  glmnet::cv.glmnet(x = xmat, y = ymat, alpha=a, nfolds=nfolds, foldid=foldid, parallel=parallel, ...)
}

#########################################################################
# Extract Optimal Alpha and Lambda Corresponding to the Lowest CV Error #
#########################################################################

extractGlmnetInfo <- function(object)
{
  # Find Lambdas
  lambda.min <- object$lambda.min
  lambda.1se <- object$lambda.1se
  
  # Determine Optimal Lambda
  which.min <- which(object$lambda == lambda.min)
  which.1se <- which(object$lambda == lambda.1se)
  
  # Summarize Optimal Lambda and Alpha
  data.frame(lambda.min = lambda.min, cv.min = object$cvm[which.min],
             lambda.1se = lambda.1se, cv.1se = object$cvm[which.1se])
}


