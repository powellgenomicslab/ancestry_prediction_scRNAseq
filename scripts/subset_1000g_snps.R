#!/usr/bin/env Rscript
.libPaths("/usr/local/lib/R/site-library")

library(data.table)
library(tidyverse)


### Set up directories ###
args <- commandArgs(trailingOnly=FALSE)
common_snps_file <- args[1]
snps_file_1000g <- args[2]
outdir <- args[3]



print("Reading in common snps file.")
common_snps <- fread(common_snps_file, sep = "_")
colnames(common_snps) <- c("chr", "bp", "a1", "a2")


print("Reading in 1000g snps file.")
snps_1000g <- fread(snps_file_1000g, sep = "_")
colnames(snps_1000g) <- c("chr", "bp", "id", "a1", "a2")


### Subset the 1000g snps by the common snps that will be used across all sites ###
subset_snps_1000g <- snps_1000g[common_snps, on = c("chr", "bp", "a1", "a2")]


snps_out_1000g <- data.table(ID = paste(subset_snps_1000g$chr, subset_snps_1000g$bp, subset_snps_1000g$id, subset_snps_1000g$a1, subset_snps_1000g$a2, sep = "_"))

fwrite(snps_out_1000g, paste0(outdir,))