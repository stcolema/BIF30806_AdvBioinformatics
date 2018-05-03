#!/usr/bin/env Rscript

# Author Stephen Coleman
# Student ID 940-309-160-050

# Initial set up of Sleuth - only required once.
# source("http://bioconductor.org/biocLite.R")
# biocLite("rhdf5")
# install.packages("devtools") 
# library("devtools")
# devtools::install_github("pachterlab/sleuth")
# biocLite("biomaRt")

library("sleuth")
library("shiny")
# Set current working directory
setwd('.')

args <- commandArgs(trailingOnly=T)

# Ensure at least one input
if(length(args) == 0){
  stop("At least one argument (input directory) must be supplied.n", call=F)
}

base_dir <- args[1]

# Want metadata of samples, similar to assignment from DESeq2 tutorial
metadata <- read.csv(args[2], header=T)
metadata <- metadata[order(metadata$run.accession),]

sample_id <- dir(file.path(base_dir))
sample_id <- grep('_quant', sample_id, value=TRUE)

entries <- c()
for(sample in metadata$run.accession){
    entry <- grep(sample, sample_id, value=T)
    entries <- c(entries, entry)
}

sample_id <- sort(entries)
kal_dirs <- sapply(entries, function(id) file.path(base_dir, id))

s2c <- dplyr::select(metadata, sample = run.accession, Tissue, SRP)
s2c <- dplyr::mutate(s2c, path = kal_dirs)

# Read the Kallisto output files, connect them with metadata, 
# and set up a linear model for analyzing the expression data.
so <- sleuth_prep(s2c, ~Tissue+SRP, extra_bootstrap_summary = TRUE)

# Next we fit the linear model and test for one of the model coefficients
# Attempt first a comparison of Tissue vs Experiment, if this fails use
# Tissue vs random
so <- tryCatch({
          so <- sleuth_fit(so, ~Tissue+SRP, 'full')
          so <- sleuth_fit(so, ~SRP, 'reduced')
      },
      error=function(cond) {
          message("Direct comparison of Tissue and Experiment failed.")
          message(cond)
          message(": Attempting alternative model.")
          so <- sleuth_fit(so, ~Tissue, 'full')
          so <- sleuth_fit(so, ~1, 'reduced')
          return(so)
}
)
# Test reduced model, full model
so <- sleuth_lrt(so, 'reduced', 'full')

sleuth_table <- sleuth_results(so, 'reduced:full', 'lrt', show_all = FALSE)
sleuth_significant <- dplyr::filter(sleuth_table, qval <= 0.05)
write.csv(sleuth_significant[1:500,], 'significant_isoforms.csv')

so <- sleuth_wt(so, which_beta="TissueLeaf") 

# Visualise the results
saveRDS(so, file = './so.rds')
sleuth_live(so, launch.browser=F)
