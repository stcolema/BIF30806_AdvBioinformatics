#!/usr/bin/env Rscript

# Author Stephen Coleman
# Student ID 940-309-160-050

# Read in arguments (currently using directory where salmon output is)
args = commandArgs(trailingOnly=TRUE)

# test if there is at least one argument: if not, return an error
if (length(args)==0) {
  stop("At least one argument must be supplied (input file).n", call.=FALSE)
}

# Various tools required
#source("http://bioconductor.org/biocLite.R")
#biocLite("devtools")    # only if devtools not yet installed
#biocLite("COMBINE-lab/wasabi")

library(wasabi)

# Use all quant files in the inputted directory to prepare for sleuth
base_dir <- file.path(args[1])
sample_id <- dir(file.path(base_dir))
sample_id <- grep('quant', sample_id, value=TRUE)
sfdirs <- sapply(sample_id, function(id) file.path(args[1], id))

# If reading in file directly uncomment the below line
#sfdirs <- file.path(args[1])

# Prepare for sleuth, hiding warnings (as file size gives warnings)
# Delete diretories that fail (as they are empty)
for(dir in sfdirs){
    tryCatch({
              prepare_fish_for_sleuth(dir)
    },
    error=function(cond) {
            unlink(dir, recursive=T)
            print(cat('Deleting ', dir, ' directory'))
            return(NA)
    }
)
}
