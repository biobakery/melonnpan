#!/usr/bin/env Rscript

################
# Load library #
################
library(melonnpan)
library(optparse)
library(GenABEL)

# Command line parameters
option_list = list(
  make_option(
    c("-i", "--metag"), # i stands for input 
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
    c("-o", "--output"), # o stands for output 
    type = "character")) 

# Print progress message
cat("Running MelonnPan-Predict using the following parameters:", "\n");
opt <- parse_args(OptionParser(option_list=option_list), positional_arguments = TRUE)
print(opt)

# Extract arguments
metag<- opt$options$metag 
weight.matrix<-opt$options$weight.matrix 
train.metag <- opt$options$train.metag 
criticalpoint <- opt$options$criticalpoint 
corr.method <- opt$options$corr.method 
output<-opt$options$output 

# Predict metabolites
DD<-melonnpan::melonnpan.predict(metag = metag, 
                             weight.matrix = weight.matrix, 
                             train.metag = train.metag,
                             criticalpoint = criticalpoint, 
                             corr.method = corr.method, 
                             output = output)
