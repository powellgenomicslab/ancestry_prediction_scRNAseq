Background
==========================

.. _freebayes: https://github.com/freebayes/freebayes
.. _data_1000g: https://www.internationalgenome.org/data/


Motivation for this Pipeline
------------------------------
Genetic ancestry is an important piece of information when assessing genetic and transcriptomic data across multiple different individuals.

Often, single cell RNA-sequencing data is produced on samples from donors whose personal information is unknown to the researchers.
For example, ancestral inforamtion is often not ascertained at sample collection and self-reporting can sometimes be misleading.
Single nucleotide polymorphism (SNP) genotyping of excess sample can be used to estimate ancestry.
However, excess sample is not always available - especially when using publicly available data.
Of course,

With this in mind, we established this pipeline to estimate genetic ancestry from scRNA-seq data by:

#. Estimating genetic information from the scRNA-seq reads (using freebayes_)
#. Aligning the samples to `1000 Genomes <data_1000g>`__ principal component (PC) space 
#. Predicting the ancestry by training a k nearest neighbors model on the `1000 Genomes <data_1000g>`__ data.


.. admonition:: Note
  :class: seealso

  We have provided instructions to run this pipeline in two ways:

  #. Through a :doc:`snakemake pipeline <RunthePipeline/index>` (suggested especially for datasets with multiple pools and individuals)

  #. :doc:`Manually <RunManually/index>` which can be used when data do not fit the assumptions in the snakemake pipeline or just to get a better idea of the steps that are implemented in the snakemake pipeline.


First Steps
-------------
First, proceed to the :doc:`Installation<Install>` section to download the required singularity image and set up the pipeline locally.




Support
-----------------
If you have any questions, suggestions or issues with any part of the Ancestry Prediction from scRNA-seq Data Pipeline, feel free to submit an `issue <https://github.com/powellgenomicslab/ancestry_prediction_scRNAseq/issues>`_ or email Drew Neavin (d.neavin @ garvan.org.au)
