#!/usr/bin/env python
import os
import pandas as pd
from glob import glob
import re
import subprocess
import numpy as np

def matchFolders(x, scrnaseq_dir, dir_list = None):
    for folder in dir_list:
        if os.path.isdir(os.path.join(scrnaseq_dir,folder)):
            if re.search(r'^' + x + "$", folder):
                return(folder)
            elif re.search(r'^' + x + '\D', folder):
                return(folder)
            elif re.search(x + "$", folder):
                return(folder)
            elif re.search(x + "\D", folder):
                return(folder)

def matchFiles(x, barcodes_dir):
    for file in os.listdir(barcodes_dir):
        if os.path.isfile(os.path.join(barcodes_dir,file)):
            if re.search(r'^' + x + "$", file):
                return(file)
            elif re.search(r'^' + x + '\D', file):
                return(file)
            elif re.search(x + "$", file):
                return(file)
            elif re.search(x + "\D", file) :
                return(file)


def get_bam_files(pool_dir):
    for dirpath, dirnames, filenames in os.walk(pool_dir):
        for filename in [f for f in filenames if f.endswith(".bam")]:
            return(os.path.join(dirpath, filename))

def get_bam_dirs(scrnaseq_filelist, pools = None):
    try:
        bam_filelist = [get_bam_files(pool) for pool in scrnaseq_filelist]
        bam_filedict = dict(zip(pools, bam_filelist))
        bamlibs = pd.Series(bam_filedict, name="Bam_Files")
        return(bamlibs)
    except Exception as error:
        print(error)
        raise SystemExit("Could not find a bam file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain '.bam' within the name.")
    


def get_scrnaseq_dirs(config):
    # Extract variables from configuration file for use within the rest of the pipeline
    input_dict = config["inputs"]
    output_dict = config["outputs"]
    ref_dict = config["refs"]

    # General variables used by the rest of the pipeline
    scrnaseq_dir = input_dict["scRNAseq_dir"]
    barcodes_dir = input_dict["barcode_annotation_dir"]

    ### Check that all the directories exist and can be accessed
    if not os.path.exists(scrnaseq_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(scrnaseq_dir))

    if not os.path.exists(barcodes_dir):
        raise Exception("Directory {} does not exist or you have not mounted a parent directory for the singularity bucket".format(barcodes_dir))

    # Read in samplesheet from the configfile specified from the command line
    samples = pd.read_csv(input_dict["metadata_file"], sep = "\t")

    # Expect first colunn to be pools
    pools = np.unique(samples.iloc[:, 0])

    # Match pools to scrna seq directories to make a list of each scRNA-seq dir
    scrna_seq_dirlist = os.listdir(scrnaseq_dir)
    try:
        scrnaseq_filelist = [os.path.join(scrnaseq_dir, matchFolders(pool, dir_list = scrna_seq_dirlist, scrnaseq_dir = scrnaseq_dir)) for pool in pools]
    except TypeError:
        print("Could not find a scRNA-seq directory for all of the pools in your pool list. Please check that they are spelled correctly and you do not have any additional pool names.")
    scrnaseq_filedict = dict(zip(pools, scrnaseq_filelist))
    scrnaseq_libs = pd.Series(scrnaseq_filedict, name="scRNAseq_Directories")


    # Match individuals to files in the barcodes directory
    try:
        annotated_barcode_filelist = [os.path.join(barcodes_dir, matchFiles(pool, barcodes_dir = barcodes_dir)) for pool in pools]
    except TypeError:
        print("Could not find an annotated barcode file for all of the pools in your metadata file in directory '" + barcodes_dir, "'. Please check that they are spelled correctly and you do not have any additional pool names.")
    annotated_barcode_filedict = dict(zip(pools, annotated_barcode_filelist))
    annotated_barcode_libs = pd.Series(annotated_barcode_filedict, name="Annotated_Barcode_Files")


    ### Get the bam files for each pool
    try:
        bam_filelist = [get_bam_files(pool) for pool in scrnaseq_filelist]
    except:
        print("Could not find a bam file in all the scRNA-seq pool directories. Please check that they exist somewhere in your pool scRNA-seq directories and contain '.bam' within the name.")
    bam_filedict = dict(zip(pools, bam_filelist))
    bamlibs = pd.Series(bam_filedict, name="Bam_Files")




    dataframe=pd.concat([scrnaseq_libs, bamlibs, annotated_barcode_libs], axis=1)
    return(dataframe)

def parsePaths(config):
    # Strip trailing slashes from paths 
    # scrnaseq_dir, individual_list_dir, output_dir
    config["inputs"]["scRNAseq_dir"] = (config["inputs"]["scRNAseq_dir"]).rstrip("/")
    config["inputs"]["barcode_annotation_dir"] = (config["inputs"]["barcode_annotation_dir"]).rstrip("/")
    config["outputs"]["outdir"] = (config["outputs"]["outdir"]).rstrip("/")
    return(config)