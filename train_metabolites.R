#!/usr/bin/env Rscript

##################
# Load libraries #
##################

library(melonnpan)
library(optparse)

###########################
# Command line parameters #
###########################

option_list = list(
  make_option(
    c("-i", "--metab"), # i stands for input metabolites 
    type = "character"),
  make_option(
    c("-g", "--metag"), # g stands for input genes or genomic features 
    type = "character"),
  make_option(
    c("-o", "--output"), # o stands for output
    type = "character"),
  make_option(
    c("-a", "--alpha"), 
    default = seq(0.05, 0.95, 0.05), 
    type = "numeric"), # a stands for alpha 
  make_option(
    c("-l", "--lambda.choice"), default="lambda.1se", # l stands for lambda
    type = "character"),
  make_option(
    c("-n", "--nfolds"), default = 10, # n stands for n-fold
    type = "integer"),
  make_option(
    c("-m", "--correction"), default = "fdr", # m stands for multiplicity correction
    type = "character"),
  make_option(
    c("-c", "--method"), default = "spearman", # c stands for correlation method
    type = "character"),
  make_option(
    c("-p", "--cores"), default = 4, # p stands for number of parallel cores 
    type = "integer"),
  make_option(
    c("-r", "--seed"), default = 1234, # r stands for random seed 
    type = "numeric"),
  make_option(
    c("-t", "--cutoff"), default = 0.3, # t stands for threshold for Cohen correlation cutoff
    type = "numeric"),
  make_option(
    c("-v", "--verbose"), default = FALSE, # v stands for verbose 
    action = "store_true"),
  make_option(
    c("-b", "--no.transform.metab"), default = FALSE, 
    action = "store_true"),
  make_option(
    c("-d", "--no.transform.metag"), default = FALSE, 
    action = "store_true"),
  make_option(
    c("-q", "--discard.poor.predictions"), default = FALSE, 
    action = "store_true"),
  make_option(
    c("-f", "--plot"), default=FALSE, # f stands for figures
    action = "store_true"),
  make_option(
    c("-s", "--outputString"), 
    default = c("MelonnPan_Training_Summary", 
                "MelonnPan_Trained_Weights",
                "MelonnPan_Trained_Metabolites"), 
    type = "character") # s stands for string
) 

##########################
# Print progress message #
##########################

cat("Running MelonnPan-Train using the following parameters:", "\n");
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
print(opt)


#########################
# Extract CLI arguments #
#########################

metab<- opt$options$metab 
metag<- opt$options$metag 
output<- opt$options$output 
alpha<- opt$options$alpha 
lambda.choice<- opt$options$lambda.choice 
nfolds<- opt$options$nfolds 
correction<- opt$options$correction 
method<- opt$options$method 
cores<- opt$options$cores 
seed<- opt$options$seed 
cutoff<- opt$options$cutoff 
verbose<- opt$options$verbose 
no.transform.metab<- opt$options$no.transform.metab 
no.transform.metag<- opt$options$no.transform.metag 
discard.poor.predictions<- opt$options$discard.poor.predictions
plot<- opt$options$plot 
outputString<- opt$options$outputString 

#####################
# Train metabolites #
#####################

DD<-melonnpan::melonnpan.train(metab = metab, 
                           metag = metag, 
                           output = output,
                           alpha = alpha, 
                           lambda.choice = lambda.choice, 
                           nfolds = nfolds,
                           correction = correction, 
                           method = method, 
                           cores = cores, 
                           seed = seed,
                           cutoff = cutoff, 
                           verbose = verbose, 
                           no.transform.metab = no.transform.metab,
                           no.transform.metag = no.transform.metag,
                           discard.poor.predictions = discard.poor.predictions,
                           plot = plot,
                           outputString = outputString)
