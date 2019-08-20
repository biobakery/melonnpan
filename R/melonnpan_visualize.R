#' Visualization of Significant Results.
#'
#' This function produces scatter plots for significantly predicted compounds.
#' @param metab Measured metabolite relative abundances (normalized to proportional data ranging from 0 to 1.0). 
#' Should have same features and samples as pred.
#' @param pred Predicted metabolite relative abundances (normalized to proportional data ranging from 0 to 1.0). 
#' Should have same features and samples as metab.
#' @param Output_file Path to the file to write the output.
#' @param cohenCorrelationCutoff Correlation threshold for significant prediction. Default is 0.3.
#' @param nrow Number of rows in batch scatter plot. Default is 2.
#' @param ncol Number of columns in batch scatter plot. Default is 2.

#' @keywords metabolite prediction, microbiome, metagenomics, elastic net, metabolomics
#' @export
melonnpan.visualize<-function(metab, 
                       pred, 
                       cohenCorrelationCutoff,
                       Output_file, 
                       ncol=2, 
                       nrow=2){
  
  ##############################################################################################
  # Calculate Correlation Coefficients and P-values between Observed and Predicted Metabolites #
  ##############################################################################################
  # This is not optimal at the moment (can be made more efficient)
  
  # Dimensionality Check
  if(all(names(metab)==names(pred))==FALSE)
    stop("Both measured and predicted tables should have the same names.")
  
  # Proportionality Check
  if(any(metab<0)||any(metab>1)|| any(pred<0)||any(pred>1))
    stop("All measurements should be normalized to proportions.")
  
  # Initialize
  scatterplot.list <- list()
  correlations <- c()
  for (m in 1:ncol(pred)) {
    corr <- cor(as.numeric(as.vector(pred[,m])), as.numeric(metab[,m]), method="spearman")
    correlations <- c(correlations,corr)
    sp <- as.data.frame(cbind(as.numeric(as.vector(pred[,m])), as.numeric(metab[,m]))) 
    colnames(sp) <- c("Predicted","Measured")
    scatterplot.list[[m]] <- ggplot(sp, aes(Predicted, Measured)) +
      ggtitle(paste(colnames(metab)[m], ": Correlation ", round(corr, 2), sep='')) + 
      geom_point(size=8, fill='purple', color="black", shape = 21, stroke = 1)  + 
      geom_smooth(method="lm", se=FALSE, color='red') +
      theme_bw() + 
      theme( panel.grid = element_blank(),  legend.position = 'none',
             axis.title.y = element_text(size = 10),
             axis.title.x = element_text(size = 10),
             axis.text.x = element_text(size = 10),
             plot.title = element_text( size = 10)) + guides(fill = FALSE) 
    rm(sp) 
  }
  which<-which(correlations>cohenCorrelationCutoff)
  scatterplot.list<-scatterplot.list[which]
  ml<-marrangeGrob(scatterplot.list, ncol=ncol, nrow=nrow,  top = NULL)
  ggsave(Output_file, ml)
  dev.off()
}
