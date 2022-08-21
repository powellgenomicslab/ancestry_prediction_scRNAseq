
Installation
==========================

This pipeline to estimate genetic ancestry from scRNA-seq data has been built in snakemake and packaged in a Singularity image that contains all the required softwares.
This should provide continuity to execution of this pipeline across different computational systems.
We're hoping that installation and pipeline execution will be relatively painless as well.


Singularity Image
--------------------
The only thing to note before you download this image is that the image is **~6.5Gb** so, depending on the internet speed, it will take **~15-30 min to download**.

To download the singularity image:

.. code-block:: bash

  wget https://www.dropbox.com/s/eqm82gv91qjwvki/ancestry_prediction_scRNAseq.sif
  wget https://www.dropbox.com/s/tubxzl3zmm3cwgk/ancestry_prediction_scRNAseq.sif.md5


Then you should check to make sure that the image downloaded completely by comparing the image md5sum to the original md5sum.
You can do that by running the following commands:

.. code-block:: bash

  md5sum ancestry_prediction_scRNAseq.sif > downloaded_ancestry_prediction_scRNAseq.sif.md5
  diff -s downloaded_ancestry_prediction_scRNAseq.sif.md5 ancestry_prediction_scRNAseq.sif.md5


If everything was downloaded correctly, that command should report:

.. code-block:: bash

  Files ancestry_prediction_scRNAseq.sif.md5 and downloaded_ancestry_prediction_scRNAseq.sif.md5 are identical

.. admonition:: Note
  :class: seealso

  Please note that the singularity image and this documentation is updated with each release. 
  This means that the most recent documentation may not be 100% compatible with the singularity image that you have.
  
  You can check the version of your singularity image to match with documentation with:

  .. code-block:: bash

    singularity inspect ancestry_prediction_scRNAseq.sif


If you run into any issues with downloading the image or any issue with running anything from this image, you can reach out to us by submitting an issue at `Github <https://github.com/powellgenomicslab/ancestry_prediction_scRNAseq/issues>`__


.. dropdown:: :octicon:`eye` Software versions - for the curious
  :color: muted



  Image build date: 24 July, 2022
 
  +----------------------------+---------------------------+-------------------------------+
  | Software Group             | Software                  | Version                       |
  +============================+===========================+===============================+
  | Supporting Softwares       | ``sinto``                 | 0.8.4                         |
  |                            +---------------------------+-------------------------------+
  |                            | ``Crossmap``              | 0.6.4                         |
  |                            +---------------------------+-------------------------------+
  |                            | ``vartrix``               | v1.1.3                        |
  |                            +---------------------------+-------------------------------+
  |                            | ``htslib``                | v1.13                         |
  |                            +---------------------------+-------------------------------+
  |                            | ``samtools``              | v1.13                         |
  |                            +---------------------------+-------------------------------+
  |                            | ``bcftools``              | v1.13                         |
  |                            +---------------------------+-------------------------------+
  |                            | ``freebayes``             | v1.3.5                        |
  +----------------------------+---------------------------+-------------------------------+
  | R Supporting Packages      | ``argparse``              | v2.1.6                        |
  | (R v4.2.1)                 +---------------------------+-------------------------------+
  |                            | ``ComplexHeatmap``        | v2.12.0                       |
  |                            +---------------------------+-------------------------------+
  |                            | ``data.table``            | v1.14.2                       |
  |                            +---------------------------+-------------------------------+
  |                            | ``vcfR``                  | v1.13.0                       |
  |                            +---------------------------+-------------------------------+
  |                            | ``tidyverse``             | v1.3.2                        |
  |                            +---------------------------+-------------------------------+
  |                            | ``cowplot``               | v1.1.1                        |
  |                            +---------------------------+-------------------------------+
  |                            |   ``colorspace``          | v2.0-3                        |
  |                            +---------------------------+-------------------------------+
  |                            |   ``ggplot2``             | v3.3.6                        |
  |                            +---------------------------+-------------------------------+
  |                            |   ``caret``               | v6.0-92                       |
  |                            +---------------------------+-------------------------------+
  |                            |   ``RColorBrewer``        | v1.1-3                        |   
  +----------------------------+---------------------------+-------------------------------+
  | Python Supporting Packages | ``argparse``              | v1.4.0                        |
  | (Python v3.6.8)            +---------------------------+-------------------------------+
  |                            | ``pysam``                 | v0.19.1                       |
  |                            +---------------------------+-------------------------------+
  |                            | ``pandas``                | v1.1.5                        |
  |                            +---------------------------+-------------------------------+
  |                            | ``scipy``                 | v1.5.4                        |
  +----------------------------+---------------------------+-------------------------------+



.. .. _common_snps:

.. Common SNP locations
.. ------------------------

.. You will need a list of common SNPs to indicate where freebayes should search for variants in the bam.
.. We have provided common SNP location bed files that can be used for calling SNPs with freebayes filtered by minor allele frequency.
.. The files contain SNPs on either hg19/GRCh37 or hg38/GRCh38 and either have 'chr' encoding or not

.. They are available to be downloaded with the links here:

.. +----------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
.. | Genome               | .. centered:: Chr Encoding   | .. centered:: vcf File                                                                                                                                            | .. centered:: md5sum File                                                                                                                                                 |
.. |                      |                              |                                                                                                                                                                   |                                                                                                                                                                           |
.. |                      |                              |                                                                                                                                                                   |                                                                                                                                                                           |
.. +======================+==============================+===================================================================================================================================================================+===========================================================================================================================================================================+
.. | GRCh37               |  .. centered:: No 'chr'      | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf <https://www.dropbox.com/s/rce452zoawd0eee/GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf>`__             | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 <https://www.dropbox.com/s/jrfb287hux6ehtg/GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5>`__             |
.. |                      +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
.. |                      | .. centered:: 'chr' encoding | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf <https://www.dropbox.com/s/qrn6df1i1wxukxn/GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf>`__ | .. centered:: `GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5 <https://www.dropbox.com/s/on47saot3d2cgij/GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5>`__ |
.. +----------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
.. | GRCh38               |  .. centered:: No 'chr'      | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf <https://www.dropbox.com/s/4nmm344g4j7pou4/GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf>`__             | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 <https://www.dropbox.com/s/izwp3l8oqwrt9dn/GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5>`__             |
.. |                      +------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+
.. |                      | .. centered:: 'chr' encoding | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf <https://www.dropbox.com/s/ycfxs407sqgoori/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf>`__ | .. centered:: `GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5 <https://www.dropbox.com/s/77opxja1cgaq994/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.vcf.md5>`__ |
.. +----------------------+------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------+---------------------------------------------------------------------------------------------------------------------------------------------------------------------------+


.. You can either just click the link for the file you want to download or you can right click > "Copy Link Address" and use wget on the command line.
.. For example, hg19/GRCh37 without 'chr' encoding:

.. .. code-block:: bash

..   wget https://www.dropbox.com/s/rce452zoawd0eee/GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf
..   wget https://www.dropbox.com/s/jrfb287hux6ehtg/GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5


.. Then check that the md5sum matches the downloaded md5sum:

.. .. code-block:: bash

..   md5sum GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf > downloaded_GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5
..   diff -s GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 downloaded_GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5


.. If everything was downloaded correctly, that command should report:

.. .. code-block:: bash

..   Files GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 and downloaded_GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.vcf.md5 are identical



Next Steps
-------------
The next section :doc:`Data Preparation<DataPreparation>` will explain the input files required for this software and their expected formats.



Support
-----------------
If you have any questions, suggestions or issues with any part of the Ancestry Prediction from scRNA-seq Data Pipeline, feel free to submit an `issue <https://github.com/powellgenomicslab/ancestry_prediction_scRNAseq/issues>`_ or email Drew Neavin (d.neavin @ garvan.org.au)
