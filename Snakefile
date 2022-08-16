#!/usr/local/envs/py36/bin python3

import os
import sys
import pandas as pd
from glob import glob
import subprocess
import shutil
import pysam
import gzip
import re
import numpy as np
from itertools import repeat

# Import custom functions
from mods import prepareArguments
from mods import read10x



### Extract variables from configuration file for use within the rest of the pipeline
config = prepareArguments.parsePaths(config)
input_dict = config["inputs"]
output_dict = config["outputs"]
ref_dict = config["refs"]
snp_dict = config["snp"]
bind_path = input_dict["bind_path"]


if (ref_dict["genome"] == "hg19"):
    fasta = ref_dict["hg19_fasta"]
elif (ref_dict["genome"] == "hg38"):
    fasta = ref_dict["hg38_fasta"]
    fasta_lift = ref_dict["hg19_fasta"]
else:
    logger.info("Your genome entry in the refs group in the yaml was not valid. It was '" + ref_dict["genome"] + "' but should be hg19 or hg38. Quitting.\n\n")




# Rule-specific arguments
freebayes_ancestry_dict = config["freebayes_ancestry"]
reference_ancestry_predictions_dict = config["reference_ancestry_predictions"]


### Check that found the correct files and report each ###
if ((ref_dict["genome"] == "hg19" and os.path.exists(fasta)) or (ref_dict["genome"] == "hg38" and os.path.exists(fasta) and os.path.exists(fasta_lift))):

    if os.path.exists(input_dict["scRNAseq_dir"]):

        if os.path.exists(input_dict["metadata_file"]):

            if os.path.exists(input_dict["singularity_image"]):

                if os.path.exists(input_dict["barcode_annotation_dir"]):

                    if not os.path.exists(output_dict["outdir"]):

                        logger.info("It appears that the output directory that you indicated does not yet exist, creating it at: {}".output_dict["outdir"])
                        os.makedirs(output_dict["outdir"])


                    # Use prepareArguments.py script to retrieve exact directories of single cell files
                    logger.info("\nReading in sample info")

                    scrnaseq_libs_df = prepareArguments.get_scrnaseq_dirs(config)
                    scrnaseq_libs_df.to_csv(os.path.join(output_dict["outdir"],'file_directories.txt'), sep = "\t", index = False)
                    samples = pd.read_csv(input_dict["metadata_file"], sep = "\t")
                    samples.columns = ["Pool", "Individual"]
                    samples['Pool'].apply(str)
                    samples['Individual'].apply(str)


                    logger.info("Done.\n\n")

                    if all([os.path.isfile(f) for f in scrnaseq_libs_df["Bam_Files"]]):


                        ### Check that chr encoding is the same between fasta reference and bam file
                        def check_if_string_in_file(file_name, first_string_to_search, string_to_search):
                            if file_name.endswith(".gz"):
                                with gzip.open(file_name, 'rt') as read_obj:
                                    lines = []
                                    match = []
                                    # Read all lines in the file one by one
                                    for line in read_obj:
                                        # For each line, check if line contains the string
                                        if re.match(first_string_to_search, line):
                                            match.append(True)
                                            if re.match(string_to_search, line):
                                                lines.append(True)
                                            break
                                    return(lines)
                            else:
                                # Open the file in read only mode
                                with open(file_name) as read_obj:
                                    lines = []
                                    # Read all lines in the file one by one
                                    for line in read_obj:
                                        # For each line, check if line contains the string
                                        if re.match(first_string_to_search, line):
                                            if re.match(string_to_search, line):
                                                lines.append(True)
                                            break
                                    return(lines)
                            if not match:
                                logger.info("Can't find the string '" + first_string_to_search + "' in your " + file_name + " file.\nPlease check this file and rerun once fixed.")

                        logger.info("Checking for chr encoding in fasta or bam files.")

                        ## Fasta check
                        fasta_chr = check_if_string_in_file(fasta, '>', '>chr')

                        if fasta_chr:
                            logger.info("   Fasta has chr encoding.")
                        else:
                            logger.info("   Fasta does not have chr encoding.")


                        ## bam check
                        bam_header = pysam.view(scrnaseq_libs_df["Bam_Files"][0], "-H")
                        if (bam_header.find("SN:chr") != -1):
                            bam_chr = True
                        else:
                            bam_chr = False

                        if bam_chr:
                            logger.info("   Bam has chr encoding.")
                        else:
                            logger.info("   Bam does not have chr encoding.")

                        if all([fasta_chr, bam_chr]) or not any([fasta_chr, bam_chr]):
                            if all([fasta_chr, bam_chr]):
                                logger.info("   Detected chr encoding in bam and fasta and will move forward assuming chr encoding for chromosomes.\n\n")
                                chr = True
                            elif not any([fasta_chr, bam_chr]):
                                logger.info("   Detected no chr encoding in bam and fasta. Will move forward assuming no chr encoding for chromosomes.\n\n")
                                chr = False


                            ### Define jobs and chain file if liftover required
                            if  ref_dict["genome"] in ["hg19", "hg38"]:
                                if ref_dict["genome"] == "hg38":
                                    if chr:
                                        chain_cross = '/opt/ancestry_prediction_scRNAseq/refs/hg38ToHg19.over.chain'
                                        bed_dir = '/opt/ancestry_prediction_scRNAseq/refs/split_beds/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding'
                                        chrs = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
                                    else:
                                        chain_cross = '/opt/ancestry_prediction_scRNAseq/refs/GRCh38_to_GRCh37.chain'
                                        bed_dir = '/opt/ancestry_prediction_scRNAseq/refs/split_beds/GRCh38_1000G_MAF0.01_GeneFiltered_NoChr'
                                        chrs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]
                                elif chr:
                                    bed = '/opt/ancestry_prediction_scRNAseq/refs/split_beds/GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding'
                                    chrs = ["chr1","chr2","chr3","chr4","chr5","chr6","chr7","chr8","chr9","chr10","chr11","chr12","chr13","chr14","chr15","chr16","chr17","chr18","chr19","chr20","chr21","chr22"]
                                else:
                                    bed = '/opt/ancestry_prediction_scRNAseq/refs/split_beds/GRCh37_1000G_MAF0.01_GeneFiltered_NoChr'
                                    chrs = [1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22]

                                # Import individual rules
                                include: "includes/freebayes_ancestry.smk"


                                #### Determine if have common snps across sites to run ####
                                post_common_snps_rules = []

                                if not (input_dict["common_snps"] in ["None", 'none', "NONE", "Null", "null", "NULL"]):

                                    logger.info("Running the remainder of the pipeline as you indicated that your snps common across all sites are listed in: " + input_dict["common_snps"])
                                    post_common_snps_rules.append(expand(output_dict["outdir"]  + "/{pool}/individual_{individual}/freebayes/pca_sex_checks_original/Ancestry_PCAs.png", zip, pool=samples.Pool, individual=samples.Individual))
                                    post_common_snps_rules.append(expand(output_dict["outdir"]  + "/{pool}/individual_{individual}/freebayes/pca_sex_checks_original/ancestry_assignments.tsv", zip, pool=samples.Pool, individual=samples.Individual))
                                    post_common_snps_rules.append(output_dict["outdir"] + "/ancestry_assignments.tsv")


                                #### Determine if need to run the reference rules and add them if necessary ####
                                reference_rules = []

                                if snp_dict["ref_snp_predict"] in ["True", "TRUE", "true", "T"]:

                                    logger.info("You indicated that you have reference SNP genotypes that you would like to use for ancestry prediction for pipeline accuracy checking.")

                                    if os.path.exists(snp_dict["vcf"]):
                                        include: "includes/reference_ancestry_predictions.smk"

                                        reference_rules.append(output_dict["outdir"] + "/reference/pca_sex_checks_original/Ancestry_PCAs.png")
                                        reference_rules.append(output_dict["outdir"] + "/reference/pca_sex_checks_original/ancestry_assignments.tsv")
                                        reference_rules.append(expand(output_dict["outdir"] + "/{pool}/souporcell/Genotype_ID_key.txt", pool=samples.Pool))
                                        reference_rules.append(output_dict["outdir"] + "/ref_sc_ancestry_prediction_comparison/assigfdccnments_probabilities_w_ref_identified.png")

                                    else:
                                        logger.info("Could not find the provided snp vcf file at: '" + snp_dict["vcf"] + "'.")

                                else:
                                    logger.info("Beginning ancestry prediction using scRNA-seq data.")

                                rule all:
                                    input:
                                        expand(output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/snps_data.tsv", zip, pool=samples.Pool, individual=samples.Individual),
                                        output_dict["outdir"] +  "/common_snps_across_pools.tsv",
                                        post_common_snps_rules,
                                        reference_rules


                            else:
                                logger.info("Your genome entry in the refs group in the yaml was not valid. It was '" + ref_dict["genome"] + "' but should be hg19 or hg38. Exiting.\n\n")

                        else:
                            logger.info("Didn't detect consistent chr chromosome encoding.\nPlease check that your bam and fasta all have the same chr encoding. Exiting.\n\n")

                    else:
                        logger.info("Could not find all the bam files from parent directory: '" + input_dict["scRNAseq_dir"] + "'.\nExiting.")

                else:
                    logger.info("Could not find the annotated barcodes file directory at: '" + input_dict["barcode_annotation_dir"] + "'.\nExiting.")

            else:
                logger.info("Could not find the singularity image at: '" + input_dict["singularity_image"] + "'.\nExiting.")

        else:
            logger.info("Could not find the provided metadata file at: '" + input_dict["metadata_file"] + "'.\nExiting.")

    else:
        logger.info("Could not find the single cell RNAseq parent directory at: '" + ref_dict["scRNAseq_dir"] + "'.\nExiting.")

else:
    logger.info("Could not find fasta file(s). Please check that you have provided the correct paths.\nExiting.")




