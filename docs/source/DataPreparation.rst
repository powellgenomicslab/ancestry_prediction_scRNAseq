Data Preparation
==========================

This section will explain the required inputs for this ancestry prediction pipeline from scRNA-seq data and the expected formats for those files.


Data Files
-----------
Here is a list of the files that you will need for this pipeline with further explanation of preparation of the files below:


.. admonition:: Required
  :class: important

  - :ref:`Bam file(s) <bam>`
    - The single cell bam files

  - :ref:`hg19 reference fasta file <hg19 fasta>`

  - :ref:`Annotated barcodes file <anno barcodes>`
    - Two-column tab separated file indicating which barcodes match to each individual in the dataset

  - :ref:`Sample metadata file <sample meta>`
    - Two-column tab separated file that contains the names of the pools and the IDs of the individual in each pool


.. admonition:: Optional

  - :ref:`hg38 reference fasta file <hg38 fasta>`
    - ONLY NEEDED IF YOUR SEQUENCE DATA WAS MAPPED TO HG38!!!

  - :ref:`Reference vcf file <ref vcf>`
    - Only needed if you have reference SNP genotypes from a microarray or called from whole genome or exome sequencing data and want to compare the ancestries predicted from the reference SNP genotypes to the ancestries predicted from the scRNA-seq data.

    





Preparation of Data Files
-----------------------------


    
.. _bam:

Bam File(s) :bdg-success-line:`Required`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Bam file(s) that contain the aligned sequences from the scRNA-seq.


.. _hg19 fasta:

hg19 Reference Fasta File :bdg-success-line:`Required`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A reference hg19 fasta that uses the same chr encoding as you bam file (i.e. chromosomes are encoded as ch1 or 1)


.. _anno barcodes:

Annotated Barcodes Files :bdg-success-line:`Required`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A tab-separated file that has the barcodes in the first column and the IDs of the individual that they are assigned to in the second column. 
This file should NOT have a header.

.. admonition:: Note
  :class: seealso

  Please make sure that your individual IDs do not start with a number as there are some softwares in the pipeline that do not handle them well.


For example:

+--------------------+--------------+
| AAACCCAAGAACTGAT-1 |      K835-8  |
+--------------------+--------------+
| AAACCCAAGAAGCCTG-1 |      K1292-4 |
+--------------------+--------------+
| AAACCCAAGCAGGTCA-1 |      K1039-4 |
+--------------------+--------------+
| AAACCCAAGCGGATCA-1 |      K962-0  |
+--------------------+--------------+
| AAACCCAAGCTGCGAA-1 |      K835-8  |
+--------------------+--------------+
| AAACCCAAGGTACTGG-1 |      K1292-4 |
+--------------------+--------------+
| AAACCCAAGTCTTCCC-1 |      K835-8  |
+--------------------+--------------+
| AAACCCACAACCGCCA-1 |      K835-8  |
+--------------------+--------------+
| AAACCCACACAGTGAG-1 |      K962-0  |
+--------------------+--------------+
| AAACCCACACCCTGAG-1 |      K835-8  |
+--------------------+--------------+
| ...                |      ...     |
+--------------------+--------------+



.. _sample meta:

Sample Metadata File :bdg-success-line:`Required`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
A tab separated file that contains two columns: the first for the Pool and the second for the 

.. admonition:: Note
  :class: seealso

  Please make sure that your individual IDs do not start with a number as there are some softwares in the pipeline that do not handle them well.


+-----------------+-------------+
| Pool            | Individual  |
+=================+=============+
| RZ731_Pool8     | K1039-4     |
+-----------------+-------------+
| RZ731_Pool8     | K1292-4     |
+-----------------+-------------+
| RZ731_Pool8     | K752-4      |
+-----------------+-------------+
| RZ731_Pool8     | K835-8      |
+-----------------+-------------+
| RZ731_Pool8     | K938-0      |
+-----------------+-------------+
| RZ731_Pool8     | K962-0      |
+-----------------+-------------+



.. _ref vcf:

Reference SNP Genotypes vcf :bdg-primary-line:`Optional`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

If you have reference SNP genotypes for the individuals in your dataset from microarray or whole exome or genome sequencing,
we have build funcionality into the pipeline to estimate ancestry based on the referene genotypes and provide comparison between the reference and scRNA-seq predicted ancestry annotations.


.. _hg38 fasta:

hg38 Reference Fasta File :bdg-primary-line:`Optional`
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

ONLY NEEDED IF YOUR SEQUENCE DATA WAS MAPPED TO HG38!!!






Support
-----------------
If you have any questions, suggestions or issues with any part of the Ancestry Prediction from scRNA-seq Data Pipeline, feel free to submit an `issue <https://github.com/powellgenomicslab/ancestry_prediction_scRNAseq/issues>`_ or email Drew Neavin (d.neavin @ garvan.org.au)
