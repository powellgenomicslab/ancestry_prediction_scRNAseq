
Preparing to Execute the Pipeline
===================================

The files required to execute the snakemake pipeline (except those listed in the :doc:`Data Preparation section <../DataPreparation>`) are packaged into the Singularity image.
Some of those files need to be transfered to the local system for pipeline execution.
The following sections help transfer these files to your local system and provide instructions on how to edit thtem


Image Setup
--------------
:bdg-success-line:`Required`

To get files from the Singularity image onto your system, simply execute:

.. code-block:: bash

  singularity run --bind <absolute_directory_path> --app setup ancestry_prediction_scRNAseq.sif <absolute_directory_path>

.. admonition:: Note
  :class: seealso

  The pipeline expects certain files pulled from the image to be in the same directory as the singularity image so you will have to rerun the setup steps if you move the image


This will transfer some of the files from the Singularity image to your local system:


.. code-block:: bash

  .
  ├── ancestry_prediction_scRNAseq.yaml
  ├── includes
  │   ├── reference_ancestry_predictions.smk
  │   └── souporcell_ancestry.smk
  ├── mods
  │   └── prepareArguments.py
  ├── Snakefile
  └── snakemake.yaml


We have included a yaml (``snakemake.yaml``) to generate a conda environment that has snakemake and accompanying software required to execute the pipeline as well as the ``ancestry_prediction_scRNAseq.yaml`` which will need to be edited based on your files and system.


Install Snakemake
-----------------------------------------
:bdg-success-line:`Required`

You will need snakemake to run the pipeline. 
We highly recommend using the conda environment that we have generated as it has all the required dependencies.
You can create this environment on your system using the ``snakemake.yaml`` that we have provided.

.. code-block:: bash

  conda env create -f snakemake.yaml -n ancestry_pred_snakemake

Then, to activate this environment:

.. code-block:: bash

  conda activate ancestry_pred_snakemake

If you would prefer to install snakemake and scipy yourself, you can follow the instructions for `installing Snakemake <https://snakemake.readthedocs.io/en/stable/getting_started/installation.html>`__ and then install scipy with ``pip install scipy`` and numpy with ``pip install numpy``.



Editing the Yaml File
--------------------------------

The last step that is needed before pipeline execution is to edit the ``ancestry_prediction_scRNAseq.yaml`` file per your system and files.
We suggest copying this file to a new filename and editing the new one.

Required User Input
^^^^^^^^^^^^^^^^^^^^^^
:bdg-success-line:`Required`

The first section of the ``ancestry_prediction_scRNAseq.yaml`` file will require user input:
 
.. code-block:: yaml

  ####################################################################################
  ##### The following arguments are for indicating file locations on your system #####
  ####################################################################################
  refs:
    genome: hg38 ## hg38 or hg19; genome the sequencing data have been aligned to
    hg19_fasta: /path/to/hg19/reference/genome.fa ## Path to the reference hg19 fasta to be used for remapping for freebayes demultiplexing steps. Ideally this would be the same reference used for original mapping but any reference on the same genome with the same 'chr' encoding will do
    hg38_fasta: /path/to/hg38/reference/genome.fa ## ONLY NEEDED IF DATA ORIGINALLY MAPPED TO HG38; Path to the reference hg38 fasta to be used for remapping for freebayes demultiplexing steps. Ideally this would be the same reference used for original mapping but any reference on the same genome with the same 'chr' encoding will do

  inputs:
    metadata_file: /path/to/samples_meta.tsv ## Sample metadata file that has two columns: 'Pool' and 'Individual. The Pool should be the exact names of the parent folders for the scRNAseq output
    singularity_image: /path/to/singularity/image.sif ### The complete path to the singularity image that has all the softwares
    bind_path: /path ## List of paths to bind to Singularity. You can specify multiple directories by adding a "," between them. Eg. ${DIRECTORY1},${DIRECTORY2}. Singularity will bind the directory that you are running from + subfolders but will not be able to find anything above unless it is in this argument
    scRNAseq_dir: /path/to/scRNAseq/parent/directory ### the parent directory that has directories for each pool and the scRNA-seq output below it
    barcode_annotation_dir: /path/to/barcodes/annotation/directory ### The directory that contains each of the barcode files with per-barcode annotation. The pool name needs to be within the file name. these should be filtered to remove doublets and contain only cells assigned to an individual 
    common_snps: None ### Leave as None for first run of the pipeline. This will be the file of SNPs common across all sites and samples. This will be generated by sending your snp list files to Drew Neavin and the garvan institute (d.neavin@garvan.org.au) to create a common list of snps.
    barcode_tag: "CB"

  outputs: 
    outdir: /path/to/parent/out/dir


Please update the contents to reflect your data and your system. 
Here is a more detailed explanation of each entry:

``refs``

- ``genome`` - This is the genome that your single cell data have been aligned to. This should be either 'hg38' or 'hg19'.

- ``hg19_fasta`` - Path to an hg19 (or GRCh37) fasta file that has the same chr encoding (`i.e.` chr1 or 1) as your aligned single cell data.

- ``hg38_fasta`` - ONLY NEEDED IF DATA ORIGINALLY MAPPED TO HG38!!! If your data was aligned to hg19 (GRCh37), you can leave this field unedited. Otherwise, provide the path to an hg38 (GRCh38) fasta file.

``inputs``

- ``metadata_file`` - Path to tab-separated file that has two columns and a header. The first column is the Pool ID. This should be the same Pool IDs used for the directories of your single cell results. The second column should have the individuals in each pool that you want processed. See example in :ref:`Data Preparation <sample meta>`.

- ``singularity_image`` - Path to singularity image that you downloaded in :doc:`../Install`.

- ``bind_path`` - Path(s) to be bound for the singularity image. Singularity by default only binds the directories and files only below where you execute the command from. Therefore, it won't be able to find any files that are elsewhere on your system. Bind as many directories as you need by separating with a comma to so all the files you need can be found.

- ``scRNAseq_dir`` - The directory that contains directories for each single cell pool below it. The pool names should match those in your metadata_file. The pipeline is built to search for bam files downstream of each pool folder. You may run in to issues if you have multiple bam files.

- ``barcode_annotation_dir`` - A directory that contains files for each of the annotated barcode files prepared as described in the :ref:`Data Preparation <anno barcodes>`.

- ``common_snps`` - As part of this project, the SNPs called from multiple sites will be combined to identfiy the common SNPs across all sites. Leave this as 'None' for the first pipeline run. After common SNPs are generated, provide the path to that file.

- ``barcode_tag`` - The tag used to indicate the barcode in your bam file. For 10x, this will be "CB" but may be different for other technologies. Update as needed.

``outputs``

- ``outdir`` - The directory where you would like all results to be output.



Reference-Based Ancestry Predictions
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-primary-line:`Optional`

We also provide the functionality to predict individual ancestries from reference SNP genotypes from microarrays or whole exome or genome sequencing data.
As part of this, the pipeline will provide comparisons of the predictions between the reference and single cell data.

To implement reference-based ancestry predictions, you will need to provide information in the second section of the ``ancestry_prediction_scRNAseq.yaml`` file:

.. code-block:: yaml

  #############################################################################################################################
  ##### The following arguments are if you have reference SNPs for these individuals for method accuracy testing purposes #####
  #############################################################################################################################
  snp:
    ref_snp_predict: False ## Set to true or false depending on if have reference SNP genotype data to be predicted for pipeline accuracy testing purposes
    vcf: /path/to/unimputed/reference/vcf.vcf ## Reference SNP genotype vcf that is UNIMPUTED 

``ref_snp_predict`` should be changed to True if you would like to predict ancestries based on reference SNP data

``vcf`` should be the complete path for a vcf containing the SNP genotypes for each individual you would like to predict ancestry for.


Additional User Inputs
^^^^^^^^^^^^^^^^^^^^^^^^^^^

Each rule has different memory and thread allocations that can be altered and set by the user. 
You may want to update this if your jobs keep failing due to limited memory allocation.
The important 

This is an example fo the top few lines of this section:

.. code-block:: yaml

  #########################################################################################################################################
  ##### The following arguments are common parameters such as memory and threads that may need to be changed depending on the dataset #####
  #########################################################################################################################################
  freebayes_ancestry:
    ### Following parameters are for bam subsetting by individual - will only be used if multi-individual sample multiplexing was used
    subset_bam_memory: 4
    subset_bam_threads: 8

    ### Following parameters are for indexing the individual subset bam
    index_memory: 4
    index_threads: 2




Support
-----------------
If you have any questions, suggestions or issues with any part of the Ancestry Prediction from scRNA-seq Data Pipeline, feel free to submit an `issue <https://github.com/powellgenomicslab/ancestry_prediction_scRNAseq/issues>`_ or email Drew Neavin (d.neavin @ garvan.org.au)
