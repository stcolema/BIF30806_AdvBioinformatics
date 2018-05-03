#!/usr/bin/env Rscript
library("optparse")

option_list = list(
  make_option(c("-f", "--folder"), type="character", default=NULL,
              help="salmon output folder name", metavar="character"),
  # make_option(c("-o", "--out"), type="character", default="out.txt",
  #             help="output file name [default= %default]", metavar="character")
);

opt_parser = OptionParser(option_list=option_list);
opt = parse_args(opt_parser);

source("http://bioconductor.org/biocLite.R")
# biocLite("devtools")    # only if devtools not yet installed
biocLite("COMBINE-lab/wasabi")

library(wasabi)

# File path and file names for salmon data
# sfdirs <- file.path("data", c("samp1", "samp2", "samp3"))
sfdirs <- file.path(opt$file)
prepare_fish_for_sleuth(sfidrs)
