Submission Examples
====================

.. _issue: https://github.com/sc-eQTLgen-consortium/WG1-pipeline-QC/issues

The submission command for snakemake is dependent on the cluster infrastructure so here are a few examples submission commands to help get you started. If you find a different submission works on your cluster and you think it could help others, feel free to send it to us so that we can add it here (either open an issue_ or email d.neavin @ garvan.org.au directly)



LSF Example
------------

Here is an example of a submission for a LSF cluster. Of course, you may have to update it based on how your cluster has been set up.

.. code-block:: bash

  nohup \
  snakemake \
  --snakefile $SCEQTL_PIPELINE_SNAKEFILE \
  --configfile $SCEQTL_PIPELINE_CONFIG \
  --rerun-incomplete \
  --jobs 20 \
  --use-singularity \
  --restart-times 2 \
  --keep-going \
  --cluster \
    "bsub \
    -W 24:00 \
    -x \
    -M 10000 \
    -e $LOG \
    -o $LOG" \
  > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &


SGE Example
-----------

This is an additional example for an SGE cluster.

.. code-block:: bash

    nohup \
      snakemake \
        --snakefile $IMPUTATION_SNAKEFILE \
        --configfile $IMPUTATION_CONFIG \
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



SLURM Examples
--------------

.. code-block:: bash

  nohup \
    snakemake \
      --snakefile $IMPUTATION_SNAKEFILE \
      --configfile $IMPUTATION_CONFIG \
      --rerun-incomplete \
      --jobs 48 \
      --use-singularity \
      --restart-times 2 \
      --keep-going \
      --cluster \
         "sbatch \
         --qos debug \
         -N 1 \
         --ntasks 1 \
         --cpus-per-task 48 \
         -o $LOG/%{rule}.out \
         --export ALL" \
       > $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &


Another SLURM example where file latency causes issues with snakemakes ability to detect if a job is completed (note the ``--latency-wait`` parameter):

.. code-block:: bash

  nohup 
    snakemake \
      --snakefile $IMPUTATION_SNAKEFILE \
      --configfile $IMPUTATION_CONFIG \
      --rerun-incomplete \
      --jobs 1 \
      --use-singularity \
      --restart-times 2 \
      --keep-going \
      --latency-wait 30 \
      --cluster \
          "sbatch \
	  --qos regular \
	  -N {threads} \
	  --mem={resources.mem_per_thread_gb}G \
	  --tmp={resources.disk_per_thread_gb}G \
	  -o $LOG/{rule}.out \
	  --export ALL \
	  --time=05:59:59" \
	> $LOG/nohup_`date +%Y-%m-%d.%H:%M:%S`.log &
