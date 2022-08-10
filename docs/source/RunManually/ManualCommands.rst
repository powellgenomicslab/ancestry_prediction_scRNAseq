Manually Execute Pipeline
=============================

This section shows each of the commands to run the ancestry prediction pipeline manually.
You can run each step of the pipeline from the Singularity image downloaded in the :doc:`../Install` section which should help standardize execution of this pipeline across different computing systems.

.. admonition:: Note
  :class: seealso

  If you have multiple pools and/or multiple individuals in each pool, you will need to update the output directories to reflect this and parallelize accross each of these combinations (following bam subsetting if you multiplexed multiple individuals in a pool).
  The snakemake pipeline is designed to automatically parallelize these jobs and monitor the state of each job so we highly recommend leveraging the pipeline if you can do so with your data.

The general steps are:

#. :ref:`Divide Bam File by Individual <divide>` :bdg-primary-line:`Optional`

#. :ref:`Index Bam File(s) <index>` :bdg-success-line:`Required`

#. :ref:`Freebayes SNP Calling <freebayes>` :bdg-success-line:`Required`

#. :ref:`Lift hg38 to hg19 <lift>` :bdg-primary-line:`Optional`

#. :ref:`Convert vcf to Plink <vcf2plink>` :bdg-success-line:`Required`

#. :ref:`Identify Common SNPs Between 1000G Data and Your Data <common_snps>` :bdg-success-line:`Required`

#. :ref:`Get the SNPs in Common Across all Pools <common_snps_pools>` :bdg-success-line:`Required`

:octicon:`stop` PAUSE - send results to Drew Neavin (d.neavin @ garvan.org.au)

#. :ref:`Subset Data for Common SNPs <subset>`  :bdg-success-line:`Required`

#. :ref:`Prune the SNPs in LD <prune>`  :bdg-success-line:`Required`

#. :ref:`Filter the 1000G and Freebayes SNPs for Pruned SNPs <filter>`  :bdg-success-line:`Required`

#. :ref:`Update Chromosome IDs <chr_update>`  :bdg-success-line:`Required`

#. :ref:`Calculate PCs for 1000G <pcs>`  :bdg-success-line:`Required`

#. :ref:`Project 1000G and Freebayes Data in PCs <pc_projection>`  :bdg-success-line:`Required`

#. :ref:`Plot PC Results <plot>` :bdg-success-line:`Required`



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

.. _divide:

1. Divide Bam File
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^
:bdg-primary-line:`Optional`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  ~15-40 minutes when using 8 threads with 4G each

If you have multiple individuals per capture, you will need to do this step but if you only have one individual in your pool, you do not have to split the bam by individuals and you can proced directly to :ref:`2. Index Bam File(s)<index>`

In preparation for this step, we set some additional parameters and create the required output directory:
The parameters that we use for the command in this step are 

.. code-block:: bash

  BAM=/path/to/bam/file.bam ### Path to bam file
  ANNO_BARCODES=/path/to/annotated/barcodes.tsv ### Path to annotated barcodes
  TAG="CB"
  N=8

  mkdir -p $OUTDIR/bams

- The ``$ANNO_BARCODES`` is the annotated barcodes file described in :doc:`../DataPreparation`

- The ``$TAG`` is the tag used in your bam file to indicate cell barcodes. In 10x captures, this is 'CB' but could be different for different technologies


To divide the bam file into a single file for each individual in the pool, simply execute:

.. code-block:: bash

  singularity exec --bind $BIND $SIF sinto filterbarcodes -b $BAM -c $ANNO_BARCODES --barcodetag $TAG --outdir $OUTDIR/bams --nproc $N



.. _index:

2. Index Bam File(s)
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min


The bam file(s) need to be indexed before SNP calling with freebayes.
Of course, if your bam is already indexed, you can skip to :ref:`Freebayes SNP Calling <freebayes>`

The ``$BAM`` will be either your original bam file (if did not subset by individual in previous step) or one of the bam files subset by each individual in the pool

.. code-block:: bash

  singularity exec --bind $BIND $SIF samtools index $BAM




.. _freebayes:

3. Freebayes SNP Calling 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  ~ 12 and 36 hours with 8 to 16 threads with 4-16G each.
  The time for this step will vary greatly depending on the number of reads per capture captured and the number of cells per individual.


Freebayes will be used to call SNP genotypes from the bam file.

You will need a list of common SNPs to indicate where freebayes should search for variants in the bam.
We have provided common SNP location bed files in the Singularity image that can be used for calling SNPs with freebayes:

- ``/opt/ancestry_prediction_scRNAseq/refs/GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.bed`` - Common SNPs filtered for 1% MAF and overlapping genes on hg19/GRCh37 reference with 'chr' chromosome encoding (`i.e.` chr1 instead of 1)

- ``/opt/ancestry_prediction_scRNAseq/refs/GRCh37_1000G_MAF0.01_GeneFiltered_NoChr.bed`` - Common SNPs filtered for 1% MAF and overlapping genes on hg19/GRCh37 reference without 'chr' chromosome encoding (`i.e.` 1 instead of chr1)

- ``/opt/ancestry_prediction_scRNAseq/refs/GRCh38_1000G_MAF0.01_GeneFiltered_ChrEncoding.bed`` - Common SNPs filtered for 1% MAF and overlapping genes on hg19/GRCh37 reference with 'chr' chromosome encoding (`i.e.` chr1 instead of 1)

- ``/opt/ancestry_prediction_scRNAseq/refs/GRCh38_1000G_MAF0.01_GeneFiltered_NoChr.bed`` - Common SNPs filtered for 1% MAF and overlapping genes on hg19/GRCh37 reference without 'chr' chromosome encoding (`i.e.` 1 instead of chr1)



Define some variables to execute CrossMap to lift the data from hg38 to hg19

.. code-block:: bash

  N=8
  TARGETS=/opt/ancestry_prediction_scRNAseq/refs/GRCh37_1000G_MAF0.01_GeneFiltered_ChrEncoding.bed ## Change this to the correct chain file for your data
  FASTA=/path/to/reference/fasta.fa



Run freebayes to identify the SNP genotyeps for the individual in the bam file:

.. code-block:: bash

  singularity exec --bind $BIND,/tmp $SIF fasta_generate_regions.py $FASTA.fai {params.regions} > $OUTDIR/regions

  export TMPDIR=/tmp
  singularity exec --bind $BIND,/tmp $SIF freebayes-parallel $OUTDIR/regionsS $N -f $FASTA -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --limit-coverage 100000 --targets $TARGETS $BAM > $OUTDIR/freebayes/freebayes.vcf






.. _lift:

4. Lift hg38 to hg19 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-primary-line:`Optional`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 10 min

If your sequence data was aligned to hg38 (GRCh38), you will need to lift it to hg19 for ancestry annotation with 1000G data.

We have provided chain files in the Singularity image that can be used for lifting the data between hg38 and hg19:

- ``/opt/ancestry_prediction_scRNAseq/refs/GRCh38_to_GRCh37.chain`` - Does not contain 'chr' enchoding (`i.e.` 1 and not chr1)

- ``/opt/ancestry_prediction_scRNAseq/refs/hg38ToHg19.over.chain`` - Does contain 'chr' enchoding (`i.e.` chr1 and not 1)


Define some variables to execute CrossMap to lift the data from hg38 to hg19

.. code-block:: bash

  CHAIN=/opt/ancestry_prediction_scRNAseq/refs/GRCh38_to_GRCh37.chain ## Change this to the correct chain file for your data
  FASTA=/path/to/reference/fasta.fa


Run CrossMap to lift the data from hg19 to hg38:

.. code-block:: bash

  singularity exec --bind $BIND $SIF CrossMap.py vcf $CHAIN $OUTDIR/freebayes/freebayes.vcf $FASTA $OUTDIR/freebayes/freebayes_hg19.vcf



.. _vcf2plink:

5. Convert vcf to Plink 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min


First, convert the vcf to the plink2 pgen files. We use ``max-alleles 2`` because downstream softwares won't be able to deal with multi-allelic sites.
The ``$VCF`` will be either ``$OUTDIR/freebayes/freebayes.vcf`` (if your sequence data was mapped to hg19/GRCh37) or ``$OUTDIR/freebayes/freebayes_hg19.vcf`` (if your sequence data was mapped to hg38/GRCh38).

.. code-block:: bash

  singularity exec --bind $BIND $SIF plink2 --vcf $VCF --make-pgen --out $OUTDIR/freebayes/freebayes --max-alleles 2


Freebayes doesn't provide an ID for each SNP that it calls but that is important for downstream SNP filtering functions so we will create a pvar file that has IDs we will make from the chromosome, basepair, allele 1 and allele 2


.. code-block:: bash

  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/freebayes.pvar $OUTDIR/freebayes/freebayes.pvar_original
  singularity exec --bind $BIND $SIF sed -i 's/^chr//g' $OUTDIR/freebayes/freebayes.pvar
  singularity exec --bind $BIND $SIF grep "#" $OUTDIR/freebayes/freebayes.pvar > $OUTDIR/freebayes/freebayes_tmp.pvar
  singularity exec --bind $BIND $SIF grep -v "#" $OUTDIR/freebayes/freebayes.pvar | awk 'BEGIN{FS=OFS="\\t"}{print $1 FS $2 FS $1 "_" $2 "_" $4 "_" $5 FS $4 FS $5 FS $6 FS $7}' >> $OUTDIR/freebayes/freebayes_tmp.pvar
  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/freebayes_tmp.pvar $OUTDIR/freebayes/freebayes.pvar



.. _common_snps:

6. Identify Common SNPs Between 1000G Data and Your Data
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min


Next, we need to subset the variants for just those that are in common between the SNP genotypes called from freebayes and those called from the 1000G data.

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  The 1000G data is located in the Singularity image as demonstrated in the below commands.


Use awk to pull SNPs that are on the same chromosomes and have the same alleles


.. code-block:: bash

  singularity exec --bind $BIND $SIF awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1,$2,$4,$5];next} ($1,$2,$4,$5) in a{print $3}' $OUTDIR/freebayes/freebayes.pvar /opt/1000G/all_phase3_filtered.pvar | sed '/^$/d' > $OUTDIR/freebayes/common_snps/snps_1000g.tsv
  singularity exec --bind $BIND $SIF awk 'BEGIN{FS=OFS="\t"}NR==FNR{a[$1,$2,$4,$5];next} ($1,$2,$4,$5) in a{print $3}' /opt/1000G/all_phase3_filtered.pvar $OUTDIR/freebayes/freebayes.pvar | sed '/^$/d' > $OUTDIR/freebayes/common_snps/snps_data.tsv




.. _common_snps_pools:

7. Get the SNPs in Common Across all Pools 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 10 min


Next, we need to identify the SNPs that are in common across all the pools (and individuals if you had multiple individuals within a given pool).
If you only have one pool, you will not need to do this step


Define some variables to get the common SNPs across all the pools

.. code-block:: bash

  META=/path/to/metadata.tsv

The metadata file is the :ref:`Sample Metadata File <sample meta>` in the :doc:`../DataPreparation` documentation.

.. code-block:: bash

  singularity exec --bind $BIND $SIF Rscript /opt/ancestry_prediction_scRNAseq/scripts/common_snps.R $META $OUTDIR


:octicon:`stop` PAUSE
^^^^^^^^^^^^^^^^^^^^^^^^^^^

.. admonition:: :octicon:`stop` PAUSE

  Send the resulting SNP file (``common_snps_across_pools.tsv``) to Drew Neavin at d.neavin @ garvan.org.au so that SNPs common across all sites can be used for ancestry annotation.
  You will need to wait until you receive the file that contains common SNPs across each site.



.. _subset:

8. Subset Data for Common SNPs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 10 min


After you have received the common SNPs file across all sites (``$COMMON_SNPS``), you can subset the 1000G and freebayes SNP data to the common SNPs.
We will use ``--rm-dup force-first`` to help deal with possible duplicate entries for the same SNP called by freebayes.

.. code-block:: bash

  ### First need to subset the 1000g snps for the SNPs common to all sites, pools and individuals ###
  singularity exec --bind $BIND $SIF Rscript /opt/ancestry_prediction_scRNAseq/scripts/subset_1000g_snps.R $COMMON_SNPS $OUTDIR/freebayes/common_snps/snps_1000g.tsv $OUTDIR

  ### Subset the freebayes-called snps for the new snps ###
  singularity exec --bind $BIND $SIF plink2 --threads {threads} --pfile {params.infile} --extract {params.snps} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out} --rm-dup force-first
  singularity exec --bind $BIND $SIF plink2 --threads {threads} --pfile {params.infile_1000g} --extract {params.snps_1000g} --make-pgen --out {params.out_1000g}



.. _prune:

9. Prune the SNPs in LD 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min


Next, filter the SNPs for those that are not in linkage disequilibrium so we have unique representation for creating PCs:

.. code-block:: bash

  singularity exec --bind $BIND $SIF plink2 --threads $N --pfile $OUTDIR/freebayes/common_snps/subset_1000g \
      --indep-pairwise 50 5 0.5 \
      --out $OUTDIR/freebayes/common_snps/subset_pruned_1000g



.. _filter:

10. Filter the 1000G and Freebayes SNPs for Pruned SNPs 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min

The only variable that needs to be defined is ``$N`` which is the number of threads you would like to use for this command:

.. code-block:: bash

  singularity exec --bind $BIND $SIF plink2 --threads $N --pfile $OUTDIR/freebayes/common_snps/subset_1000g --extract $OUTDIR/freebayes/common_snps/subset_pruned_1000g.prune.out --rm-dup 'force-first' --make-pgen --out $OUTDIR/freebayes/common_snps/subset_pruned_1000g


  ### If have comments in the freebayes pvar file, need to transfer them
  if [[ $(grep "##" $OUTDIR/freebayes/common_snps/subset_data.pvar | wc -l) > 0 ]]
  then
      singularity exec --bind $BIND $SIF grep "##" $OUTDIR/freebayes/common_snps/subset_data.pvar > $OUTDIR/freebayes/common_snps/subset_pruned_data_1000g_key.txt
  fi

  singularity exec --bind $BIND $SIF awk -F"\\t" 'BEGIN{OFS=FS = "\\t"} NR==FNR{a[$1 FS $2 FS $4 FS $5] = $0; next} {ind = $1 FS $2 FS $4 FS $5} ind in a {print a[ind], $3}' $OUTDIR/freebayes/common_snps/subset_pruned_1000g.pvar $OUTDIR/freebayes/common_snps/subset_data.pvar | \
      singularity exec --bind $BIND $SIF grep -v "##" >> $OUTDIR/freebayes/common_snps/subset_pruned_data_1000g_key.txt
  singularity exec --bind $BIND $SIF grep -v "##" $OUTDIR/freebayes/common_snps/subset_pruned_data_1000g_key.txt | \
      singularity exec --bind $BIND $SIF awk 'BEGIN{FS=OFS="\t"}{print $NF}' > $OUTDIR/freebayes/common_snps/subset_data.prune.out
  singularity exec --bind $BIND $SIF plink2 --threads $N --pfile {params.infile} --extract $OUTDIR/freebayes/common_snps/subset_data.prune.out --rm-dup 'force-first' --make-pgen 'psam-cols='fid,parents,sex,phenos --out $OUTDIR/freebayes/common_snps/subset_pruned_data
  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/common_snps/subset_pruned_data.pvar $OUTDIR/freebayes/common_snps/subset_pruned_data_original.pvar
  singularity exec --bind $BIND $SIF grep -v "#" $OUTDIR/freebayes/common_snps/subset_pruned_data_original.pvar | \
      singularity exec --bind $BIND $SIF awk 'BEGIN{FS=OFS="\t"}{print($3)}' > $OUTDIR/freebayes/common_snps/SNPs2keep.txt
  singularity exec --bind $BIND $SIF grep "#CHROM" $OUTDIR/freebayes/common_snps/subset_pruned_data_1000g_key.txt > $OUTDIR/freebayes/common_snps/subset_pruned_data.pvar
  singularity exec --bind $BIND $SIF grep -Ff $OUTDIR/freebayes/common_snps/SNPs2keep.txt $OUTDIR/freebayes/common_snps/subset_pruned_data_1000g_key.txt >> $OUTDIR/freebayes/common_snps/subset_pruned_data.pvar
  singularity exec --bind $BIND $SIF awk 'BEGIN{FS=OFS="\t"}NF{NF-=1};1' < $OUTDIR/freebayes/common_snps/subset_pruned_data.pvar > $OUTDIR/freebayes/common_snps/subset_pruned_data_temp.pvar
  singularity exec --bind $BIND $SIF grep "##" $OUTDIR/freebayes/common_snps/subset_pruned_1000g.pvar > $OUTDIR/freebayes/common_snps/subset_pruned_data.pvar
  singularity exec --bind $BIND $SIF sed -i "/^$/d" $OUTDIR/freebayes/common_snps/subset_pruned_data_temp.pvar
  singularity exec --bind $BIND $SIF cat $OUTDIR/freebayes/common_snps/subset_pruned_data_temp.pvar >> $OUTDIR/freebayes/common_snps/subset_pruned_data.pvar
  singularity exec --bind $BIND $SIF plink2 --threads $N --pfile $OUTDIR/freebayes/common_snps/subset_pruned_1000g --rm-dup 'force-first' --make-bed --out $OUTDIR/freebayes/common_snps/subset_pruned_1000g 




.. _chr_update:

11. Update Chromosome IDs
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min

  
Plink has some requirements for the chromosome IDs for SNPs on the X, Y and MT chromosomes - they need to be updated to numeric.
Update the chromosome IDs to ensure compatibility with plink.

.. code-block:: bash

  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/common_snps/final_subset_pruned_data.bed $OUTDIR/freebayes/common_snps/split/
  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/common_snps/final_subset_pruned_data.bim $OUTDIR/freebayes/common_snps/split/
  singularity exec --bind $BIND $SIF sed -i 's/^X/23/g' $OUTDIR/freebayes/common_snps/split/final_subset_pruned_data.bim
  singularity exec --bind $BIND $SIF sed -i 's/^Y/24/g' $OUTDIR/freebayes/common_snps/split/final_subset_pruned_data.bim
  singularity exec --bind $BIND $SIF sed -i 's/^XY/25/g' $OUTDIR/freebayes/common_snps/split/final_subset_pruned_data.bim
  singularity exec --bind $BIND $SIF sed -i 's/^MT/26/g' $OUTDIR/freebayes/common_snps/split/final_subset_pruned_data.bim
  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/common_snps/final_subset_pruned_data.fam $OUTDIR/freebayes/common_snps/split/
  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/common_snps/subset_pruned_1000g.bed $OUTDIR/freebayes/common_snps/split
  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/common_snps/subset_pruned_1000g.bim $OUTDIR/freebayes/common_snps/split
  singularity exec --bind $BIND $SIF cp $OUTDIR/freebayes/common_snps/subset_pruned_1000g.fam $OUTDIR/freebayes/common_snps/split
  singularity exec --bind $BIND $SIF sed -i 's/^X/23/g' $OUTDIR/freebayes/common_snps/split/subset_pruned_1000g.bim
  singularity exec --bind $BIND $SIF sed -i 's/^Y/24/g' $OUTDIR/freebayes/common_snps/split/subset_pruned_1000g.bim
  singularity exec --bind $BIND $SIF sed -i 's/^XY/25/g' $OUTDIR/freebayes/common_snps/split/subset_pruned_1000g.bim
  singularity exec --bind $BIND $SIF sed -i 's/^MT/26/g' $OUTDIR/freebayes/common_snps/split/subset_pruned_1000g.bim



.. _pcs:

12. Calculate PCs for 1000G 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min


Next, calculate the principal components using the 1000G SNPs (alrleady filtered for the same SNPs as your dataset).
Again, the only additional variable that you need to define is the number of threads (``$N``) you want to use to calculate the PCs:

.. code-block:: bash

  singularity exec --bind $BIND $SIF plink2 --threads $N --pfile $OUTDIR/freebayes/common_snps/subset_pruned_1000g \
      --freq counts \
      --pca allele-wts \
      --out $OUTDIR/freebayes/pca_projection/subset_pruned_1000g_pcs



.. _pc_projection:

13. Project 1000G and Freebayes Data in PCs 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min


Next, project the 1000G and freebayes SNPs in the PCs calculated in the last step.
Again, the only additional variable that you need to define is the number of threads (``$N``) you want to use to calculate the PCs:

.. code-block:: bash

  export OMP_NUM_THREADS=$N
  singularity exec --bind $BIND $SIF plink2 --threads $N --pfile $OUTDIR/freebayes/common_snps/final_subset_pruned_data \
      --read-freq $OUTDIR/freebayes/pca_projection/subset_pruned_1000g_pcs.acount \
      --score $OUTDIR/freebayes/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
              variance-standardize \
      --score-col-nums 6-15 \
      --out $OUTDIR/freebayes/pca_projection/final_subset_pruned_data_pcs
  singularity exec --bind $BIND $SIF plink2 --threads $N --pfile $OUTDIR/freebayes/common_snps/subset_pruned_1000g \
      --read-freq $OUTDIR/freebayes/pca_projection/subset_pruned_1000g_pcs.acount \
      --score $OUTDIR/freebayes/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele 2 5 header-read no-mean-imputation \
              variance-standardize \
      --score-col-nums 6-15 \
      --out $OUTDIR/freebayes/pca_projection/subset_pruned_1000g_pcs_projected



.. _plot:

14. Plot PC Results 
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

:bdg-success-line:`Required`

.. admonition:: :octicon:`stopwatch` Expected Timing
  :class: seealso

  < 5 min


Lastly, we can predict sample SNP-based ancestry and produce some figures for visualization using our wrapper script:

.. code-block:: bash

  singularity exec --bind $BIND $SIF echo $OUTDIR/freebayes/pca_sex_checks_original/ > $OUTDIR/freebayes/pca_sex_checks_original/variables.tsv
  singularity exec --bind $BIND $SIF echo $OUTDIR/freebayes/pca_projection/final_subset_pruned_data_pcs.sscore >> $OUTDIR/freebayes/pca_sex_checks_original/variables.tsv
  singularity exec --bind $BIND $SIF echo $OUTDIR/freebayes/pca_projection/subset_pruned_1000g_pcs_projected.sscore >> $OUTDIR/freebayes/pca_sex_checks_original/variables.tsv
  singularity exec --bind $BIND $SIF echo $OUTDIR/freebayes/common_snps/subset_1000g.psam >> $OUTDIR/freebayes/pca_sex_checks_original/variables.tsv
  singularity exec --bind $BIND $SIF Rscript /opt/ancestry_prediction_scRNAseq/scripts/PCA_Projection_Plotting_original.R $OUTDIR/freebayes/pca_sex_checks_original/variables.tsv

