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
    c("-i", "--metag"), # i stands for input 
    type = "character"),
  make_option(
    c("-o", "--output"), # o stands for output 
    type = "character"),
  make_option(
    c("-w", "--weight.matrix"), default=NULL, # w stands for weight 
    type = "character"),
  make_option(
    c("-r", "--train.metag"), default=NULL, # r stands for RTSI
    type = "character"),
  make_option(
    c("-s", "--criticalpoint"), default = 0.9793, # s stands for significance threshold
    type = "numeric"),
  make_option(
    c("-c", "--corr.method"), default = "pearson", # c stands for correlation
    type = "character"),
  make_option(
    c("-b", "--no.transform.metab"), default = FALSE, 
    action = "store_true"),
  make_option(
    c("-d", "--no.transform.metag"), default = FALSE, 
    action = "store_true")) 

##########################
# Print progress message #
##########################

cat("Running MelonnPan-Predict using the following parameters:", "\n");
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
print(opt)

#####################
# Extract arguments #
#####################

metag<- opt$options$metag 
output<-opt$options$output 
weight.matrix<-opt$options$weight.matrix 
train.metag <- opt$options$train.metag 
criticalpoint <- opt$options$criticalpoint 
corr.method <- opt$options$corr.method 
no.transform.metab<- opt$options$no.transform.metab 
no.transform.metag<- opt$options$no.transform.metag 

######################################
# Predict metabolites in new samples #
######################################

DD<-melonnpan::melonnpan.predict(metag = metag, 
                                 output = output,
                                 weight.matrix = weight.matrix, 
                                 train.metag = train.metag,
                                 criticalpoint = criticalpoint, 
                                 corr.method = corr.method,
                                 no.transform.metab = no.transform.metab,
                                 no.transform.metag = no.transform.metag)
                   
