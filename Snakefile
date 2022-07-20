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
souporcell_ancestry_dict = config["souporcell_ancestry"]
reference_ancestry_predictions_dict = config["reference_ancestry_predictions"]


### Check that found the correct files and report each ###
if ((ref_dict["genome"] == "hg19" and os.path.exists(fasta)) or (ref_dict["genome"] == "hg38" and os.path.exists(fasta) and os.path.exists(fasta_lift))):

    if os.path.exists(ref_dict["vcf"]):

        if os.path.exists(input_dict["metadata_file"]):

            if os.path.exists(input_dict["singularity_image"]):

                if os.path.exists(output_dict["outdir"]):


                    # Use prepareArguments.py script to retrieve exact directories of single cell files
                    logger.info("\nReading in sample info")

                    scrnaseq_libs_df = prepareArguments.get_scrnaseq_dirs(config)
                    scrnaseq_libs_df.to_csv(os.path.join(output_dict["outdir"],'file_directories.txt'), sep = "\t", index = False)
                    samples = pd.read_csv(input_dict["metadata_file"], sep = "\t")
                    samples.columns = ["Pool", "N"]

                    samples = pd.read_csv(input_dict["metadata_file"], sep = "\t")

                    logger.info("Done.\n\n")

                    if all([os.path.isfile(f) for f in scrnaseq_libs_df["Bam_Files"]]):


                        ### Check that chr encoding is the same between fasta reference, vcf and bam file
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

                        logger.info("Checking for chr encoding in fasta, vcf and bam files.")

                        ## Fasta check
                        fasta_chr = check_if_string_in_file(fasta, '>', '>chr')

                        if fasta_chr:
                            logger.info("   Fasta has chr encoding.")
                        else:
                            logger.info("   Fasta does not have chr encoding.")


                        ## vcf check
                        vcf_chr = check_if_string_in_file(ref_dict["vcf"], '##contig=<ID', 'chr')

                        if vcf_chr:
                            logger.info("   Vcf has chr encoding.")
                        else:
                            logger.info("   Vcf does not have chr encoding.")


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

                        if all([fasta_chr, vcf_chr, bam_chr]) or not any([fasta_chr, vcf_chr, bam_chr]):
                            if all([fasta_chr, vcf_chr, bam_chr]):
                                logger.info("   Detected chr encoding in bam, vcf and fasta and will move forward assuming chr encoding for chromosomes.\n\n")
                                chr = True
                            elif not any([fasta_chr, vcf_chr, bam_chr]):
                                logger.info("   Detected no chr encoding in bam, vcf and fasta. Will move forward assuming no chr encoding for chromosomes.\n\n")
                                chr = False


                            ### Define jobs and chain file if liftover required
                            if  ref_dict["genome"] in ["hg19", "hg38"]:
                                if ref_dict["genome"] == "hg38":
                                    if chr:
                                        chain_cross = '/opt/ancestry_prediction_scRNAseq/refs/hg38ToHg19.over.chain'
                                    else:
                                        chain_cross = '/opt/ancestry_prediction_scRNAseq/refs/GRCh38_to_GRCh37.chain'


                                # Import individual rules
                                include: "includes/souporcell_ancestry.smk"



                                #### Determine if need to run the reference rules and add them if necessary ####
                                reference_rules = []

                                if snp_dict["ref_snp_predict"] == True:

                                    logger.info("You indicated that you have reference SNP genotypes that you would like to use for ancestry prediction for pipeline accuracy checking.")

                                    if os.path.exists(snp_dict["vcf"]):
                                        include: "includes/reference_ancestry_predictions.smk"

                                        reference_rules.append(output_dict["outdir"] + "/reference/pca_sex_checks_original/Ancestry_PCAs.png")
                                        reference_rules.append(output_dict["outdir"] + "/reference/pca_sex_checks_original/ancestry_assignments.tsv")
                                        reference_rules.append(expand(output_dict["outdir"] + "/{pool}/souporcell/Genotype_ID_key.txt", pool=samples.Pool))
                                        reference_rules.append(output_dict["outdir"] + "/ref_sc_ancestry_prediction_comparison/assignments_probabilities_w_ref_identified.png")

                                    else:
                                        logger.info("Could not find the provided snp vcf file at: '" + snp_dict["vcf"] + "'.")

                                else:
                                    logger.info("Beginning ancestry prediction using scRNA-seq data.")

                                rule all:
                                    input:
                                        expand(output_dict["outdir"] + "/{pool}/souporcell/clusters.tsv", pool=samples.Pool),
                                        expand(output_dict["outdir"] + "/{pool}/souporcell/pca_sex_checks_original/Ancestry_PCAs.png", pool=samples.Pool),
                                        expand(output_dict["outdir"] + "/{pool}/souporcell/pca_sex_checks_original/ancestry_assignments.tsv", pool=samples.Pool),
                                        reference_rules


                            else:
                                logger.info("Your genome entry in the refs group in the yaml was not valid. It was '" + ref_dict["genome"] + "' but should be hg19 or hg38. Quitting.\n\n")

                        else:
                            logger.info("Didn't detect consistent chr chromosome encoding.\nPlease check that your bam, vcf and fasta all have the same chr encoding. Quitting.\n\n")

                    else:
                        logger.info("Could not find all the bam files from parent directory: '" + input_dict["scRNAseq_dir"] + "'.\nExiting.")

                else:
                    logger.info("Could not find the output directory at: '" + output_dict["outdir"] + "'.\nExiting.")

            else:
                logger.info("Could not find the singularity image at: '" + input_dict["singularity_image"] + "'.\nExiting.")

        else:
            logger.info("Could not find the provided metadata file at: '" + input_dict["metadata_file"] + "'.\nExiting.")

    else:
        logger.info("Could not find the provided reference vcf file at: '" + ref_dict["vcf"] + "'.\nExiting.")

else:
    logger.info("Could not find fasta file(s). Please check that you have provided the correct paths.\nExiting.")




