Annotate Ancestry of Samples Using Reference Genotypes
==========================================================

We have provided instructions here to annotate your samples that have genotype data on your samples with either SNP genotype data or whole genome or exome sequencing.

This will allow users that have additional genotype data for their samples to compare the predictions from the sc sequencihng to the reference genotype data.


To do this, you will need your data will need to be on the GRCh37 (hg19) reference - sorry, we haven't provided instructions or resources for predicting ancestry for data directly on GRCh38 (hg38).
So if your data is on GRCh38 (hg38), you will have to first lift your datat to hg19 (GRCh37).


The general steps are:

#. :ref:`Convert Data to Plink pgen Format <convert_pgen>` :bdg-primary-line:`Optional`

#. :ref:`Identify Common SNPs Between 1000G Data and Your Data <common_snps_ref>` :bdg-success-line:`Required`

#. :ref:`Prune the SNPs in LD <prune_ref>`  :bdg-success-line:`Required`

#. :ref:`Filter the 1000G and Freebayes SNPs for Pruned SNPs <filter_ref>` :bdg-success-line:`Required`

#. :ref:`Update Chromosome IDs <chr_update_ref>`  :bdg-success-line:`Required`

#. :ref:`Calculate PCs for 1000G <pcs_ref>`  :bdg-success-line:`Required`

#. :ref:`Project 1000G and Freebayes Data in PCs <pc_projection_ref>`  :bdg-success-line:`Required`

#. :ref:`Plot PC Results <plot_ref>` :bdg-success-line:`Required`

#. :ref:`Compare Reference to single cell-predicted Ancestry <compare_ref>` :bdg-primary-line:`Optional`



Setup
---------

There are a few variables that will be used throughout the example commands which are best to define in a file that you can easily run at each step or source for execution of each step.

- ``$BIND`` - the path(s) on your system to bind to singularity when executing commands. By default, Singularity binds just the directory (and downstream dierectory and files) from your current working directory when executing the command. However, if you have some files elsewhere on your system, you can provide a parent directory to the singulairty command to indicate which additional directories to bind. Multiple directories can be included and separated by a comma (`i.e.` $DIR1,$DIR2)

- ``$SIF`` - the path to the singularity image downloaded in the :doc:`../Install` section

- ``$OUTDIR`` - the path to the output directory where all results will be written and saved

.. code-block:: bash

  BIND = /bind/path 
  SIF = /path/to/singularity/image/ancestry_prediction_scRNAseq.sif
  OUTDIR = /path/to/base/outdir





Steps
---------

.. _convert_pgen:

1. Convert Data to Plink pgen Format
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-primary-line:`Optional`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 10 minutes but may be longer with large datasets

This step will convert your file to a Plink pgen format which is required for input into the downstream steps.
If your data is already in the Plink pgen format, you can skip this step and move directly to :ref:`Identify Common SNPs Between 1000G Data and Your Data <common_snps>`.

In preparation for this step, we set some additional parameters and create the required output directory.
The parameters that we use for the command in this step are:


.. tab-set::

  .. tab-item:: VCF File
    :sync: key1

    .. code-block:: bash
        
      VCF=/path/to/vcf.vcf

    - This is pretty self-explanatory - this is the path to your sample genotype file.


  .. tab-item:: Plink BED Files
    :sync: key2

    .. code-block:: bash
        
      BED_BASE=/path/to/plink/bed_basename

    - The ``$BED_BASE`` is the  your sample genotype Plink bedfiles basename. For example if your files are named ``data.bed``, ``data.fam`` and ``data.bim`` your ``bed_basename`` would be ``data``.


.. tab-set::

  .. tab-item:: VCF File
    :sync: key1

    .. code-block:: bash
              
      singularity exec --bind $BIND $SIF plink2 --vcf $VCF --make-pgen --out $OUTDIR/data_name

    .. admonition:: Note
      :class: seealso

      This command will generate three files using the out basename provided: ``$OUTDIR/data_name.pgen``, ``$OUTDIR/data_name.pvar``, ``$OUTDIR/data_name.psam``



  .. tab-item:: Plink BED Files
    :sync: key2

    .. code-block:: bash

      singularity exec --bind $BIND $SIF plink2 --bfile $BED_BASE --make-pgen --out $OUTDIR/data_name

    .. admonition:: Note
      :class: seealso

      This command will generate three files using the out basename provided: ``$OUTDIR/data_name.pgen``, ``$OUTDIR/data_name.pvar``, ``$OUTDIR/data_name.psam``



.. _common_snps_ref:

2. Identify Common SNPs Between 1000G Data and Your Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-success-line:`Required`

This step will identify the genetic variants in common between your data and the 1000G reference and then subset your data and 1000G reference data for the genetic variants in common between the two datasets.
The 1000G reference data is provided within the Singularity image so you will not have to download it separately.
The paths to the 1000G reference data are the paths within the Singularity image.


.. code-block:: bash

  singularity exec --bind $BIND $SIF awk 'NR==FNR{a[$1,$2,$4,$5];next} ($1,$2,$4,$5) in a{print $3}' $OUTDIR/data_name.pvar /opt/1000G/all_phase3_filtered.pvar > $OUTDIR/common_snps/snps_1000g.tsv
  singularity exec --bind $BIND $SIF awk 'NR==FNR{a[$1,$2,$4,$5];next} ($1,$2,$4,$5) in a{print $3}' /opt/1000G/all_phase3_filtered.pvar $OUTDIR/data_name.pvar > $OUTDIR/common_snps/snps_data.tsv

  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile $OUTDIR/data_name --extract $OUTDIR/common_snps/snps_data.tsv --make-pgen 'psam-cols='fid,parents,sex,phenos --out $OUTDIR/common_snps/subset_data
  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile /opt/1000G/all_phase3_filtered --extract $OUTDIR/common_snps/snps_1000g.tsv --make-pgen --out $OUTDIR/common_snps/subset_1000g



.. _prune_ref:

3. Prune the SNPs in LD
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-success-line:`Required`

This step will prune the genotype data so that the genetic variants that are in low linkage disequilibrium.

.. code-block:: bash

  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile /opt/1000G/all_phase3_filtered \
      --indep-pairwise 50 5 0.5 \
      --out $OUTDIR/common_snps/subset_pruned_1000g



.. _filter_ref:

4. Filter the 1000G and Freebayes SNPs for Pruned SNPs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-success-line:`Required`

This step will prune the genotype data so that the genetic variants that are in low linkage disequilibrium (uniquely representing genetic variation across the genome).

.. code-block:: bash

  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile /opt/1000G/all_phase3_filtered --extract $OUTDIR/common_snps/subset_pruned_1000g.prune.out --make-pgen --out $OUTDIR/common_snps/subset_pruned_1000g

  if [[ $(grep "##" $OUTDIR/common_snps/subset_1000g.pvar | wc -l) > 0 ]]
  then
      singularity exec --bind $BIND $SIF grep "##" $OUTDIR/common_snps/subset_1000g.pvar > $OUTDIR/common_snps/subset_pruned_data_1000g_key.txt
  fi

  singularity exec --bind $BIND $SIF awk -F"\\t" 'BEGIN{OFS=FS = "\\t"} NR==FNR{a[$1 FS $2 FS $4 FS $5] = $0; next} {ind = $1 FS $2 FS $4 FS $5} ind in a {print a[ind], $3}' $OUTDIR/common_snps/subset_pruned_1000g.pvar $OUTDIR/common_snps/subset_pruned_data.pvar | singularity exec --bind $BIND $SIF grep -v "##" >> $OUTDIR/common_snps/subset_pruned_data_1000g_key.txt
  singularity exec --bind $BIND $SIF grep -v "##" $OUTDIR/common_snps/subset_pruned_data_1000g_key.txt | singularity exec --bind $BIND $SIF awk 'BEGIN{FS=OFS="\t"}{print $NF}' > $OUTDIR/common_snps/subset_data.prune.out
  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile $OUTDIR/common_snps/subset_data --extract $OUTDIR/common_snps/subset_data.prune.out --make-pgen 'psam-cols='fid,parents,sex,phenos --out $OUTDIR/common_snps/subset_pruned_data
  singularity exec --bind $BIND $SIF cp $OUTDIR/common_snps/subset_pruned_data.pvar $OUTDIR/common_snps/subset_pruned_data_original.pvar
  singularity exec --bind $BIND $SIF grep -v "#" $OUTDIR/common_snps/subset_pruned_data_original.pvar | singularity exec --bind $BIND $SIF awk 'BEGIN{FS=OFS="\t"}{print($3)}' > $OUTDIR/common_snps/SNPs2keep.txt
  singularity exec --bind $BIND $SIF grep "#CHROM" $OUTDIR/common_snps/subset_pruned_data_1000g_key.txt > $OUTDIR/common_snps/subset_pruned_data.pvar
  singularity exec --bind $BIND $SIF grep -Ff $OUTDIR/common_snps/SNPs2keep.txt $OUTDIR/common_snps/subset_pruned_data_1000g_key.txt >> $OUTDIR/common_snps/subset_pruned_data.pvar
  singularity exec --bind $BIND $SIF awk 'BEGIN{FS=OFS="\t"}NF{NF-=1};1' < $OUTDIR/common_snps/subset_pruned_data.pvar > $OUTDIR/common_snps/subset_pruned_data_temp.pvar
  singularity exec --bind $BIND $SIF grep "##" $OUTDIR/common_snps/subset_pruned_data_temp.pvar > $OUTDIR/common_snps/subset_pruned_data.pvar
  singularity exec --bind $BIND $SIF cat $OUTDIR/common_snps/subset_pruned_data_temp.pvar >> $OUTDIR/common_snps/subset_pruned_data.pvar
  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile $OUTDIR/common_snps/subset_pruned_1000g --make-bed --out $OUTDIR/common_snps/subset_pruned_1000g

  ### This is a contingency to remove duplicated snps from both 
  singularity exec --bind $BIND $SIF plink2 --rm-dup 'force-first' -threads 2 --pfile $OUTDIR/common_snps/subset_data --make-pgen 'psam-cols='fid,parents,sex,phenos --out $OUTDIR/common_snps/final_subset_pruned_data
  singularity exec --bind $BIND $SIF plink2 --rm-dup 'force-first' -threads 2 --pfile $OUTDIR/common_snps/subset_data --make-bed --out $OUTDIR/common_snps/subset_pruned_data


.. _chr_update_ref:

5. Update Chromosome IDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-success-line:`Required`

Plink requires specific chromosome encodings for chromosomes X, Y, and mitochondria so we will update them to be sure.

.. code-block:: bash

  singularity exec --bind $BIND $SIF cp $OUTDIR/common_snps/final_subset_pruned_data.bed $OUTDIR/common_snps/split/
  singularity exec --bind $BIND $SIF cp $OUTDIR/common_snps/final_subset_pruned_data.bim $OUTDIR/common_snps/split/
  singularity exec --bind $BIND $SIF sed -i 's/^X/23/g' $OUTDIR/common_snps/split/final_subset_pruned_data.bim
  singularity exec --bind $BIND $SIF sed -i 's/^Y/24/g' $OUTDIR/common_snps/split/final_subset_pruned_data.bim
  singularity exec --bind $BIND $SIF sed -i 's/^XY/25/g' $OUTDIR/common_snps/split/final_subset_pruned_data.bim
  singularity exec --bind $BIND $SIF sed -i 's/^MT/26/g' $OUTDIR/common_snps/split/final_subset_pruned_data.bim
  singularity exec --bind $BIND $SIF cp $OUTDIR/common_snps/final_subset_pruned_data.fam $OUTDIR/common_snps/split/
  singularity exec --bind $BIND $SIF cp $OUTDIR/common_snps/subset_pruned_1000g.bed $OUTDIR/common_snps/split/
  singularity exec --bind $BIND $SIF cp $OUTDIR/common_snps/subset_pruned_1000g.bim $OUTDIR/common_snps/split/
  singularity exec --bind $BIND $SIF cp $OUTDIR/common_snps/subset_pruned_1000g.fam $OUTDIR/common_snps/split/
  singularity exec --bind $BIND $SIF sed -i 's/^X/23/g' $OUTDIR/common_snps/split/subset_pruned_1000g.bim
  singularity exec --bind $BIND $SIF sed -i 's/^Y/24/g' $OUTDIR/common_snps/split/subset_pruned_1000g.bim
  singularity exec --bind $BIND $SIF sed -i 's/^XY/25/g' $OUTDIR/common_snps/split/subset_pruned_1000g.bim
  singularity exec --bind $BIND $SIF sed -i 's/^MT/26/g' $OUTDIR/common_snps/split/subset_pruned_1000g.bim


.. _pcs_ref:

6. Calculate PCs for 1000G
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-success-line:`Required`

Next, we will calculate principal components using the 1000G reference data.

.. code-block:: bash

  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile $OUTDIR/common_snps/subset_pruned_1000g \
      --freq counts \
      --pca allele-wts \
      --out $OUTDIR/pca_projection/subset_pruned_1000g_pcs



.. _pc_projection_ref:

7. Project 1000G and Freebayes Data in PCs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-success-line:`Required`

Next, let's project the data (both your data and the 1000G reference data) into the principal component space that you calculated in the last step.

.. code-block:: bash

  export OMP_NUM_THREADS=2

  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile $OUTDIR/common_snps/final_subset_pruned_data \
      --read-freq $OUTDIR/pca_projection/subset_pruned_1000g_pcs.acount \
      --score $OUTDIR/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
              variance-standardize \
      --score-col-nums 6-15 \
      --out $OUTIR/pca_projection/final_subset_pruned_data_pcs

  singularity exec --bind $BIND $SIF plink2 --threads 2 --pfile {params.infile_1000g} \
      --read-freq $OUTDIR/pca_projection/subset_pruned_1000g_pcs.acount \
      --score $OUTDIR/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
              variance-standardize \
      --score-col-nums 6-15 \
      --out $OUTDIR/pca_projection/subset_pruned_1000g_pcs_projected


.. _plot_ref:

8. Plot PC Results
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-success-line:`Required`

Fnally, let's plot the results and predict the ancestry of the samples in your dataset.

.. code-block:: bash

  singularity exec --bind $BIND $SIF echo $OUTDIR/pca_sex_checks_original > $OUTDIR/pca_sex_checks_original/variables.tsv
  singularity exec --bind $BIND $SIF echo $OUTDIR/pca_projection/final_subset_pruned_data_pcs.sscore >> $OUTDIR/pca_sex_checks_original/variables.tsv
  singularity exec --bind $BIND $SIF echo $OUTDIR/pca_projection/subset_pruned_1000g_pcs_projected.sscore >> $OUTDIR/pca_sex_checks_original/variables.tsv
  singularity exec --bind $BIND $SIF echo $OUTDIR/common_snps/subset_1000g.psam >> $OUTDIR/pca_sex_checks_original/variables.tsv
  singularity exec --bind $BIND $SIF Rscript /opt/ancestry_prediction_scRNAseq/scripts/PCA_Projection_Plotting_original.R $OUTDIR/pca_sex_checks_original/variables.tsv


.. _compare_ref:

9. Compare Reference to single cell-predicted Ancestry
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-primary-line:`Optional`

The last step is to compare the predicted ancestries in your genotyped data (either microarray or whole genome or exome sequencing data) to the predicted ancestries from your single cell data. 
We assume that the ID you use for the sample is the same in bpoth the single cell and genotyped data.

You will need to provide a metadata file for each of the Pools and samples in each pool as input for this step.
The file is the same :ref:`Sample metadata file <sample meta>`.

In preparation for this step, we set some additional parameters and create the required output directory.
The parameters that we use for the command in this step are:

.. code-block:: bash

  META=/path/to/meta_file.tsv

We assume that the results for each of the samples that you have single cell sequencing are located in the same base directory in the following format:

.. code-block:: bash

  .
  ├── Pool1
  │   ├── individual_1
  │   ├── ...
  │   └── individual_n
  │   ...
  └── Poolm
      ├── individual_1
      ├── ...
      └── individual_n

.. code-block:: bash

  singularity exec --bind $BIND $SIF echo $OUTDIR > $OUTDIR/ref_sc_ancestry_prediction_comparison/variables.tsv
  singularity exec --bind $BIND $SIF echo $OUTDIR/pca_sex_checks_original/ancestry_assignments.tsv >> $OUTDIR/ref_sc_ancestry_prediction_comparison/variables.tsv
  singularity exec --bind $BIND $SIF echo $OUTDIR >> $OUTDIR/ref_sc_ancestry_prediction_comparison/variables.tsv
  singularity exec --bind $BIND $SIF echo $META >> $OUTDIR/ref_sc_ancestry_prediction_comparison/variables.tsv
  singularity exec --bind $BIND $SIF Rscript /opt/ancestry_prediction_scRNAseq/scripts/compare_ref_seq_snp_ancestries.R $OUTDIR/ref_sc_ancestry_prediction_comparison/variables.tsv




Results
----------
After running the final step, you should have the following results directories.

We've highlighted the key results files (``Ancestry_PCAs.png`` and ``ancestry_assignments.tsv``):

.. code-block:: bash
  :emphasize-lines: 50,51

  .
  ├── common_snps
  │   ├── final_subset_pruned_data.bed
  │   ├── final_subset_pruned_data.bim
  │   ├── final_subset_pruned_data.fam
  │   ├── final_subset_pruned_data.log
  │   ├── final_subset_pruned_data.pgen
  │   ├── final_subset_pruned_data.psam
  │   ├── final_subset_pruned_data.pvar
  │   ├── snps_1000g.tsv
  │   ├── SNPs2keep.txt
  │   ├── snps_data.tsv
  │   ├── subset_1000g.log
  │   ├── subset_1000g.pgen
  │   ├── subset_1000g.psam
  │   ├── subset_1000g.pvar
  │   ├── subset_data.log
  │   ├── subset_data.pgen
  │   ├── subset_data.prune.out
  │   ├── subset_data.psam
  │   ├── subset_data.pvar
  │   ├── subset_pruned_1000g.bed
  │   ├── subset_pruned_1000g.bim
  │   ├── subset_pruned_1000g.fam
  │   ├── subset_pruned_1000g.log
  │   ├── subset_pruned_1000g.pgen
  │   ├── subset_pruned_1000g.popu
  │   ├── subset_pruned_1000g.prune.in
  │   ├── subset_pruned_1000g.prune.out
  │   ├── subset_pruned_1000g.psam
  │   ├── subset_pruned_1000g.pvar
  │   ├── subset_pruned_data_1000g_key.txt
  │   ├── subset_pruned_data.log
  │   ├── subset_pruned_data_original.pvar
  │   ├── subset_pruned_data.pgen
  │   ├── subset_pruned_data.psam
  │   ├── subset_pruned_data.pvar
  │   └── subset_pruned_data_temp.pvar
  ├── pca_projection
  │   ├── final_subset_pruned_data_pcs.log
  │   ├── final_subset_pruned_data_pcs.sscore
  │   ├── subset_pruned_1000g_pcs.acount
  │   ├── subset_pruned_1000g_pcs.eigenval
  │   ├── subset_pruned_1000g_pcs.eigenvec
  │   ├── subset_pruned_1000g_pcs.eigenvec.allele
  │   ├── subset_pruned_1000g_pcs.log
  │   ├── subset_pruned_1000g_pcs_projected.log
  │   └── subset_pruned_1000g_pcs_projected.sscore
  ├── pca_sex_checks_original
  │   ├── ancestry_assignments.tsv
  │   ├── Ancestry_PCAs.png
  │   └── variables.tsv
  ├── reference.log
  ├── reference.pgen
  ├── reference.psam
  └── reference.pvar



- The ``Ancestry_PCAs.png`` figure shows the 1000G individual locations in PC space compared to the individuals in each pool. For example:

  .. figure:: ../_figures/ref_Ancestry_PCAs.png
    :align: center
    :figwidth: 700px
    

- The ``ref_ancestry_assignments.tsv`` file has the annotations and probabilities for each pool. For example:

+---------+---------+-----------------+---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+------------+------------+---------+---------+---------+---------+---------+---------------------+------------------+
| FID     | IID     | PC1             | PC2           | PC3             | PC4             | PC5             | PC6             | PC7             | PC8             | PC9        | PC10       | AFR     | AMR     | EAS     | EUR     | SAS     | combined_assignment | Final_Assignment |
+=========+=========+=================+===============+=================+=================+=================+=================+=================+=================+============+============+=========+=========+=========+=========+=========+=====================+==================+
| 0       | 1       | 0.137           | -0.108        | -0.025          | -0.042          | 0.032           | -0.042          | 0.001           | -0.021          | -0.086     | -0.019     | 0       | 0       | 1       | 0       | 0       | EAS                 | EAS              |
+---------+---------+-----------------+---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+------------+------------+---------+---------+---------+---------+---------+---------------------+------------------+
| ...     | ...     |  ...            | ...           | ...             | ...             | ...             | ...             | ...             | ...             | ...        | ...        | ...     | ...     | ...     | ...     | ...     | ...                 | ...              |
+---------+---------+-----------------+---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+------------+------------+---------+---------+---------+---------+---------+---------------------+------------------+



