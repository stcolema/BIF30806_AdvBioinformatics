#!/usr/bin/env Rscript
# Tutorial from SciLifeLab
# https://scilifelab.github.io/courses/rnaseq/labs/kallisto

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

# Function for connecting ENSEMBL transcript names to 
# common gene names. This will turn out to be useful 
# at the end, when we look at the dynamic visualization
# of the results.
tx2gene <- function(){
  mart <- biomaRt::useMart(biomart = "ensembl", dataset = "hsapiens_gene_ensembl")
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id", "ensembl_gene_id",
                                       "external_gene_name"), mart = mart)
  t2g <- dplyr::rename(t2g, target_id = ensembl_transcript_id,
                       ens_gene = ensembl_gene_id, ext_gene = external_gene_name)
  return(t2g)
}

t2g <- tx2gene()

base_dir <- "C:/Users/steph/Desktop/Bioinformatics/BIF30806 - Advanced Bioinformatics/Project/data/SciLife_Kallisto_Sleuth_tut"
samples <- paste0("sample", 1:12)
kal_dirs <- sapply(samples, function(id) file.path(base_dir, id))

# Metadata of samples, similar to assignment from DESeq2 tutorial
s2c <- data.frame(path=kal_dirs, sample=samples, 
                  timepoint = rep(c("ctrl", "t2h", "t6h", "t24h"),
                                  each=3
                                  ), 
                  stringsAsFactors=FALSE
                  )

# Again, if there were other experimental factors involved, 
# these could have been modelled here as well

# Read the Kallisto output files, connect them with metadata, 
# and set up a linear model for analyzing the expression data.
so <- sleuth_prep(s2c, ~timepoint, target_mapping = t2g)

# Next we fit the linear model and test for one of the model coefficients. In this case we test the 24h time point versus the control.
so <- sleuth_fit(so)
so <- sleuth_wt(so, which_beta="timepointt24h") 


# Visualise the results
sleuth_live(so)
