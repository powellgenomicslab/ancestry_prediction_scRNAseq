Ancestry Prediction Pipeline Execution
============================================

.. _issue: https://github.com/powellgenomicslab/ancestry_prediction_scRNAseq/issues
.. _Snakemake: https://snakemake.readthedocs.io/en/stable/getting_started/installation.html

If you have any questions or issues, feel free to open an issue_ or directly email Drew Neavin (d.neavin @ garvan.org.au)


This section will provide commands to run the pipeline and expected results.

A detailed explanation of each step in the pipeline is provided in the :doc:`../RunManually/ManualCommands` documents.



Pipeline Execution Phases
-------------------------------
Phase 1
^^^^^^^^^^^^
The first phase will identify the SNPs from each individual in each pool in your dataset.
Then, we ask that you send that SNP list to Drew Neavin (d.neavin @ garvan.org.au) and she will identify the variants in common across all sites.
She will send the list of SNPs common to all sites to each user.

Phase 2
^^^^^^^^^^^^
The second phase will use the list of SNPs common to all sites to predict ancestries for each individual in the dataset.

The pipeline has been set up to run to the 'pause' point and then to execute the remainder of the rules once the SNPs common to all sites has been provided.

Optional Reference Ancestry Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
In addition, we have provided additional functionality to predict SNP-based ancestry from reference SNP genotypes (from microarrays or whole exome or genome sequencing).
This also includes an automatic comparison between the reference and single cell-predicted ancestries for accuracy testing purposes.
This is not dependent on the common SNP list across sites so can be executed at any time.
See the :doc:`PreparingYaml` documentation for more information on setting up your ``ancestry_prediction_scRNAseq.yaml`` file for this functionality.



Preparation
--------------

Before running Snakemake_, let's define some variables that will make it easier and cleaner to run.

Let's define two variables, ``$SNAKEFILE_ANC`` (the location of the ``Snakefile``) and ``$CONFIG_ANC``  which is (the location of the edited ``ancestry_prediction_scRNAseq.yaml``).
These files were copied from the singularity image when you ran the setup command:

.. code-block:: bash
  :emphasize-lines: 2,8

  .
  ├── ancestry_prediction_scRNAseq.yaml
  ├── includes
  │   ├── reference_ancestry_predictions.smk
  │   └── souporcell_ancestry.smk
  ├── mods
  │   └── prepareArguments.py
  ├── Snakefile
  └── snakemake.yaml


Change the path based on where the files is on your system:

.. code-block:: bash

  SNAKEFILE_ANC=/path/to/ancestry_prediction_scRNAseq/Snakefile
  CONFIG_ANC=/path/to/edited/ancestry_prediction_scRNAseq.yaml

Finally, let's define a location that we would like to keep the cluster log outputs and make the directory.
Change the path based on where you want log files to be written.

.. code-block:: bash

  LOG=/path/to/cluster/log_dir
  mkdir -p $LOG


.. important:: 

  If you log out of the cluster and log back in to run more jobs, you will have to redefine each of those variables. 
  We recommend keeping the commands in a file that can easily be used to define each variable or putting them in a file that you can ``source`` before running Snakemake_ each time.


Running the Pipeline
-------------------------------------

Phase 1: SNP Calling
^^^^^^^^^^^^^^^^^^^^^^^^^^^^
Now we have all the pieces that we need to run the pipeline.
This Snakemake_ pipeline is built to run all the SNP genotype imputation pre-processing steps that are necessary.
As discussed previously, the pipeline will be run in two phases.
The first phase will identify the genetic variants in each sample.
If you have reference SNP genotypes (*i.e.* microarray SNP genotype data or whole genome or exome data), some of those jobs will be run here too


#. First, let's do a "dry run" to identify what jobs will be run (remember to activate you snakemake environment before running: ``conda activate ancestry_pred_snakemake``):


   .. code-block:: bash

    snakemake \
        --snakefile $SNAKEFILE_ANC \
        --configfile $CONFIG_ANC \
        --dryrun \
        --cores 1 \
        --quiet

   The result should show you all the jobs that snakemake will run:

   .. tab-set::
  
    .. tab-item:: Without Reference SNP Genotypes
      :sync: key1

      .. code-block:: bash

        Job counts:
              count   jobs
              1       all
              1       common_snps_across_pools
              132     freebayes
              6       freebayes_common_snps
              6       freebayes_merge
              6       freebayes_update_vcf
              6       freebayes_vcf2plink
              6       index
              1       subset_bam
              165

    .. tab-item:: With Reference SNP Genotypes
        :sync: key2

        .. code-block:: bash

          Job counts:
              count   jobs
              1       all
              1       common_snps_across_pools
              132     freebayes
              1       freebayes_combine_results
              6       freebayes_common_snps
              6       freebayes_final_pruning
              6       freebayes_merge
              6       freebayes_pca_1000g
              6       freebayes_pca_project
              6       freebayes_pca_projection_assign_original
              6       freebayes_prune_1000g
              6       freebayes_update_vcf
              6       freebayes_vcf2plink
              6       index
              1       reference_common_snps
              1       reference_final_pruning
              1       reference_freebayes_comparison
              1       reference_pca_1000g
              1       reference_pca_project
              1       reference_pca_projection_assign_original
              1       reference_prune_1000g
              1       reference_vcf2plink
              1       subset_bam
              6       subset_common_snps
              210

   .. admonition:: Note
    :class: seealso

    The number of rules to be run will depend on the number of samples and pools that you have.
    The number of ``subset_bam`` rules should reflect the number of pools you have.
    The number of all other rules should be the number of samples you have.



#. Next we can check how each of these jobs relates to one another:

   .. code-block:: bash

    snakemake \
        --snakefile $SNAKEFILE_ANC \
        --configfile $CONFIG_ANC \
        --dag | \
        dot -Tsvg \
            > dag1.svg

   .. tab-set::
  
    .. tab-item:: Without Reference SNP Genotypes
      :sync: key1

      The resulting image will be saved to your current directory. In this case, we are using just one pool with 6 individuals for illustration purposes but this figure will change depending on the number of pools and individuals in your dataset. There's quite a lot in this figure so if you would like to see it you can view it `here <https://user-images.githubusercontent.com/44268007/185299587-f224146c-d3ed-4f03-bf97-a4daa022ca26.svg>`__.

    .. tab-item:: With Reference SNP Genotypes
      :sync: key2

      The resulting image will be saved to your current directory. In this case, we are using just one pool with 6 individuals for illustration purposes but this figure will change depending on the number of pools and individuals in your dataset. There's quite a lot in this figure so if you would like to see it you can view it `here <https://user-images.githubusercontent.com/44268007/185781307-94dcdf14-33af-4b14-9b2f-35458dc9eeb3.svg>`__.



#. Next, let's run those jobs:

   .. admonition:: Important
    :class: caution

    You will likely need to change the cluster command dependent on your job submission platform.
    This example is the job submission command for an SGE cluster. Some other submission examples for SLURM, LSF and SGE clusters are available in :doc:`SubmissionExamples` documentation.


   .. code-block:: bash

    nohup \
      snakemake \
        --snakefile $SNAKEFILE_ANC \
        --configfile $CONFIG_ANC \
        --rerun-incomplete \
        --jobs 20 \
        --use-singularity \
        --restart-times 2 \
        --keep-going \
        --cluster \
            "qsub -S /bin/bash \
            -q short.q \
            -r yes \
            -pe smp {threads} \
            -l tmp_requested={resources.disk_per_thread_gb}G \
            -l mem_requested={resources.mem_per_thread_gb}G \
            -e $LOG \
            -o $LOG \
            -j y \
            -V" \
      > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &


.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  ~12-48 hours to run depending on the number of cells per individual and the coverage of SNPs




:octicon:`stop` PAUSE
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. admonition:: :octicon:`stop` PAUSE

  Send the resulting SNP file (``common_snps_across_pools.tsv``) which should be in your base output directory to Drew Neavin at d.neavin @ garvan.org.au so that SNPs common across all sites can be used for ancestry annotation.
  You will need to wait until you receive the file that contains common SNPs across each site.




Phase 2: Ancestry Prediction
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Preparation
+++++++++++++++++++++

After you have received the list of SNP genotypes that were identified at each site from Drew, you will have to add the location of this file to your ``ancestry_prediction_scRNAseq`` in ``common_snps:`` (highlighted below):

.. code-block:: bash
  :emphasize-lines: 15

  ####################################################################################
  ##### The following arguments are for indicating file locations on your system #####
  ####################################################################################
  refs:
    genome: hg38 ## hg38 or hg19; genome the sequencing data have been aligned to
    hg19_fasta: /path/to/hg19/reference/genome.fa ## Path to the reference hg19 fasta to be used for remapping for freebayes demultiplexing steps. Ideally this would be the same reference used for original mapping but any reference on the same genome with the same 'chr' encoding will do
    hg38_fasta: /path/to/hg38/reference/genome.fa ## ONLY NEEDED IF DATA ORIGINALLY MAPPED TO HG38; Path to the reference hg38 fasta to be used for remapping for freebayes demultiplexing steps. Ideally this would be the same reference used for original mapping but any reference on the same genome with the same 'chr' encoding will do

  inputs:
    metadata_file: /path/to/samples_meta.tsv ## Sample metadata file that has two columns: 'Pool' and 'N'. The Pool should be the exact names of the parent folders for the scRNAseq output
    singularity_image: /path/to/singularity/image.sif ### The complete path to the singularity image that has all the softwares
    bind_path: /path ## List of paths to bind to Singularity. You can specify multiple directories by adding a "," between them. Eg. ${DIRECTORY1},${DIRECTORY2}. Singularity will bind the directory that you are running from + subfolders but will not be able to find anything above unless it is in this argument
    scRNAseq_dir: /path/to/scRNAseq/parent/directory ### the parent directory that has directories for each pool and the scRNA-seq output below it
    barcode_annotation_dir: /path/to/barcodes/annotation/directory ### The directory that contains each of the barcode files with per-barcode annotation. The pool name needs to be within the file name. these should be filtered to remove doublets and contain only cells assigned to an individual 
    common_snps: /path/to/common_snps_across_sites.tsv ### Leave as None for first run of the pipeline. This will be the file of SNPs common across all sites and samples. This will be generated by sending your snp list files to Drew Neavin and the garvan institute (d.neavin@garvan.org.au) to create a common list of snps.
    barcode_tag: "CB"

  outputs: 
    outdir: /path/to/parent/out/dir


Execution
+++++++++++++++++++++

Now that we have provided the path to the SNP genotypes that will be used for ancestrys predictions, we can move on to execute the file steps of the pipeline:

#. Let's first do another dry run to see what steps will be run.

   .. code-block:: bash

    snakemake \
      --snakefile $SNAKEFILE_ANC \
      --configfile $CONFIG_ANC \
      --dryrun \
      --cores 1 \
      --reason


   - The result should show you all the jobs that snakemake will run:

   .. tab-set::
  
    .. tab-item:: Without Reference SNP Genotypes
      :sync: key1

      .. code-block:: bash

        Job counts:
          count	jobs
          1	all
          1	freebayes_combine_results
          6	freebayes_final_pruning
          6	freebayes_pca_1000g
          6	freebayes_pca_project
          6	freebayes_pca_projection_assign_original
          6	freebayes_prune_1000g
          6	subset_common_snps
          38

    .. tab-item:: With Reference SNP Genotypes
      :sync: key2

      .. code-block:: bash

        Job counts:
          count	jobs
          1	all
          1	freebayes_combine_results
          6	freebayes_final_pruning
          6	freebayes_pca_1000g
          6	freebayes_pca_project
          6	freebayes_pca_projection_assign_original
          6	freebayes_prune_1000g
          1	reference_freebayes_comparison
          6	subset_common_snps
          39

   .. admonition:: Note
    :class: seealso

    The number of rules to be run will depend on the number of samples and pools that you have.
    The number for each rule should be the number of samples that you have.
    For this example we have oine pool that has 6 total samples.




#. Let's also take a look at how those new jobs fit in with the steps that we already ran:

   .. code-block:: bash

    snakemake \
        --snakefile $SNAKEFILE_ANC \
        --configfile $CONFIG_ANC \
        --dag | \
        dot -Tsvg \
            > dag2.svg

   .. tab-set::
  
    .. tab-item:: Without Reference SNP Genotypes
      :sync: key1

      The resulting image will show jobs that are completed in dashed lines and those that still need to be run in solid lines. This will be saved to your current directory.
      The resulting saved image will show jobs that are completed in dashed lines and those that still need to be run in solid lines. 
      In this case, we are using just one pool with 6 individuals for illustration purposes but this figure will change depending on the number of pools and individuals in your dataset.
      There's quite a lot in this figure so if you would like to see it you can view it `here <https://user-images.githubusercontent.com/44268007/185772200-03ad3b5f-91a1-4bef-bf8f-9773f9f519b6.svg>`__.

    .. tab-item:: With Reference SNP Genotypes
      :sync: key2

      The resulting image will show jobs that are completed in dashed lines and those that still need to be run in solid lines. This will be saved to your current directory.
      The resulting saved image will show jobs that are completed in dashed lines and those that still need to be run in solid lines. 
      In this case, we are using just one pool with 6 individuals for illustration purposes but this figure will change depending on the number of pools and individuals in your dataset.
      There's quite a lot in this figure so if you would like to see it you can view it `here <https://user-images.githubusercontent.com/44268007/185792786-c236a793-fb9e-4d7e-8363-a57f36f0d922.svg>`__.




#. Next, let's run those new jobs:

   .. admonition:: Note
    :class: seealso

    Remember that you may need to change the cluster command dependent on your job submission platform.
    This example is the job submission command for an SGE cluster.

   .. code-block:: bash

    nohup \
      snakemake \
        --snakefile $SNAKEFILE_ANC \
        --configfile $CONFIG_ANC \
        --rerun-incomplete \
        --jobs 20 \
        --use-singularity \
        --restart-times 2 \
        --keep-going \
        --cluster \
            "qsub -S /bin/bash \
            -q short.q \
            -r yes \
            -pe smp {threads} \
            -l tmp_requested={resources.disk_per_thread_gb}G \
            -l mem_requested={resources.mem_per_thread_gb}G \
            -e $LOG \
            -o $LOG \
            -j y \
            -V" \
      > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &



Results
^^^^^^^^^^^^^

Now you have run the complete pipeline and all of your results should be in the output directory that you indicated.
An explanation of the results are in the :doc:`Results` documentation.