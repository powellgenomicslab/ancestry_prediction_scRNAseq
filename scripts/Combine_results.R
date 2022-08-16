#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")

##### Read in Libraries #####
library(tidyverse)
library(data.table)



args <- commandArgs(trailingOnly = TRUE)
outdir <- args[1]
samples_file <- args[2]




print("Reading in samples dataframe.")
samples <- fread(samples_file, sep = "\t")
colnames(samples) <- c("Pool", "Individual")


##### Read in the ancestry prediction dataframes #####
dt_list <- list()

for (row in 1:nrow(samples)){
    dt_list[[paste0(samples[row,"Pool"], "_", samples[row,"Individual"])]] <- fread(paste0(outdir, "/", samples[row,"Pool"], "/individual_", samples[row,"Individual"], "/freebayes/pca_sex_checks_original/ancestry_assignments.tsv"), sep = "\t")
}



##### Account for possibility that there may only be one individual in one pool at one location
if (nrow(samples) > 1 ){

        dt <- do.call(rbind, dt_list)

} else {

    dt <- dt_list[[1]]

}


fwrite(dt, paste0(outdir, "/ancestry_assignments.tsv"), sep = "\t")