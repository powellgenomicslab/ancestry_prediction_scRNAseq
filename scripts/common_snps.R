#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")

library(data.table)
library(tidyverse)


### Set up directories ###
args <- commandArgs(trailingOnly=FALSE)
samples_file <- args[1]
outdir <- args[2]


print("Reading in samples dataframe.")
samples <- fread(samples_file, sep = "\t")
colnames(samples) <- c("Pool", "Individual")



##### Read in the snp dataframes #####
snps_list <- list()

for (row in 1:nrow(samples)){
    snps_list[[paste0(samples[row,"Pool"], "_", samples[row,"Individual"]]] <- fread(paste0(outdir, "/", samples[row,"Pool"], "/individual_", samples[row,"Individual"], "/freebayes/common_snps/snps_data.tsv"))
}


snps_df <- snps_list %>% reduce(left_join, by = "ID")

fwrite(snps_df, paste0(outdir, "/common_snps_across_pools.tsv"))