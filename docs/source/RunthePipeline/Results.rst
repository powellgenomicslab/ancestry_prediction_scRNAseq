Results!
=========
After running those jobs, you should be done! 


Overview of the Results
--------------------------
The following results should be in your output direcotry.
The first tab provides a truncated results tree to help with visualization but you can see the complete tree from one pool with 6 samples in the second tab (Complete Tree).

.. tabs::

  .. tab:: Truncated Tree

    .. code-block:: bash

      .
      ├── ancestry_assignments.tsv
      ├── common_snps_across_pools.tsv
      ├── file_directories.txt
      ├── Pool1
      │   ├── bams
      │   │   ├── 1.bam
      │   │   ├── 1.bam.bai
      │   │   ├── 2.bam
      │   │   ├── 2.bam.bai
      │   │   ├── 3.bam
      │   │   ├── 3.bam.bai
      │   │   ├── 4.bam
      │   │   ├── 4.bam.bai
      │   │   ├── 5.bam
      │   │   ├── 5.bam.bai
      │   │   ├── 6.bam
      │   │   ├── 6.bam.bai
      │   │   └── subset_bam.done
      │   ├── individual_1
      │   │   ├── common_snps
      │   │   │   ├── final_subset_pruned_data.bed
      │   │   │   ├── final_subset_pruned_data.bim
      │   │   │   ├── final_subset_pruned_data.fam
      │   │   │   ├── final_subset_pruned_data.log
      │   │   │   ├── final_subset_pruned_data.pgen
      │   │   │   ├── final_subset_pruned_data.psam
      │   │   │   ├── final_subset_pruned_data.pvar
      │   │   │   ├── snps_1000g.tsv
      │   │   │   ├── SNPs2keep.txt
      │   │   │   ├── snps_data.tsv
      │   │   │   ├── subset_1000g.log
      │   │   │   ├── subset_1000g.pgen
      │   │   │   ├── subset_1000g.psam
      │   │   │   ├── subset_1000g.pvar
      │   │   │   ├── subset_data.log
      │   │   │   ├── subset_data.pgen
      │   │   │   ├── subset_data.prune.out
      │   │   │   ├── subset_data.psam
      │   │   │   ├── subset_data.pvar
      │   │   │   ├── subset_pruned_1000g.bed
      │   │   │   ├── subset_pruned_1000g.bim
      │   │   │   ├── subset_pruned_1000g.fam
      │   │   │   ├── subset_pruned_1000g.log
      │   │   │   ├── subset_pruned_1000g.pgen
      │   │   │   ├── subset_pruned_1000g.popu
      │   │   │   ├── subset_pruned_1000g.prune.in
      │   │   │   ├── subset_pruned_1000g.prune.out
      │   │   │   ├── subset_pruned_1000g.psam
      │   │   │   ├── subset_pruned_1000g.pvar
      │   │   │   ├── subset_pruned_data_1000g_key.txt
      │   │   │   ├── subset_pruned_data.log
      │   │   │   ├── subset_pruned_data_original.pvar
      │   │   │   ├── subset_pruned_data.pgen
      │   │   │   ├── subset_pruned_data.psam
      │   │   │   ├── subset_pruned_data.pvar
      │   │   │   └── subset_pruned_data_temp.pvar
      │   │   ├── freebayes_hg19.vcf
      │   │   ├── freebayes_hg19.vcf.unmap
      │   │   ├── freebayes.log
      │   │   ├── freebayes.pgen
      │   │   ├── freebayes.psam
      │   │   ├── freebayes.pvar
      │   │   ├── freebayes.pvar_original
      │   │   ├── freebayes_tmp.pvar
      │   │   ├── freebayes.vcf
      │   │   ├── pca_projection
      │   │   │   ├── final_subset_pruned_data_pcs.log
      │   │   │   ├── final_subset_pruned_data_pcs.sscore
      │   │   │   ├── subset_pruned_1000g_pcs.acount
      │   │   │   ├── subset_pruned_1000g_pcs.eigenval
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec.allele
      │   │   │   ├── subset_pruned_1000g_pcs.log
      │   │   │   ├── subset_pruned_1000g_pcs_projected.log
      │   │   │   └── subset_pruned_1000g_pcs_projected.sscore
      │   │   └── pca_sex_checks_original
      │   │       ├── ancestry_assignments.tsv
      │   │       ├── Ancestry_PCAs.png
      │   │       └── variables.tsv
      │   ├── individual_2
      │   │   ...
      │   ├── individual_3
      │   │   ...
      │   ├── individual_4
      │   │   ...
      │   ├── individual_5
      │   │   ├── common_snps
      │   │   ...
      │   └── individual_6
      │   └── ...
      └── snps_1000g_common_across_sites.tsv



  .. tab:: Complete Tree

    .. code-block:: bash

      .
      ├── ancestry_assignments.tsv
      ├── common_snps_across_pools.tsv
      ├── file_directories.txt
      ├── Pool1
      │   ├── bams
      │   │   ├── 1.bam
      │   │   ├── 1.bam.bai
      │   │   ├── 2.bam
      │   │   ├── 2.bam.bai
      │   │   ├── 3.bam
      │   │   ├── 3.bam.bai
      │   │   ├── 4.bam
      │   │   ├── 4.bam.bai
      │   │   ├── 5.bam
      │   │   ├── 5.bam.bai
      │   │   ├── 6.bam
      │   │   ├── 6.bam.bai
      │   │   └── subset_bam.done
      │   ├── individual_1
      │   │   ├── common_snps
      │   │   │   ├── final_subset_pruned_data.bed
      │   │   │   ├── final_subset_pruned_data.bim
      │   │   │   ├── final_subset_pruned_data.fam
      │   │   │   ├── final_subset_pruned_data.log
      │   │   │   ├── final_subset_pruned_data.pgen
      │   │   │   ├── final_subset_pruned_data.psam
      │   │   │   ├── final_subset_pruned_data.pvar
      │   │   │   ├── snps_1000g.tsv
      │   │   │   ├── SNPs2keep.txt
      │   │   │   ├── snps_data.tsv
      │   │   │   ├── subset_1000g.log
      │   │   │   ├── subset_1000g.pgen
      │   │   │   ├── subset_1000g.psam
      │   │   │   ├── subset_1000g.pvar
      │   │   │   ├── subset_data.log
      │   │   │   ├── subset_data.pgen
      │   │   │   ├── subset_data.prune.out
      │   │   │   ├── subset_data.psam
      │   │   │   ├── subset_data.pvar
      │   │   │   ├── subset_pruned_1000g.bed
      │   │   │   ├── subset_pruned_1000g.bim
      │   │   │   ├── subset_pruned_1000g.fam
      │   │   │   ├── subset_pruned_1000g.log
      │   │   │   ├── subset_pruned_1000g.pgen
      │   │   │   ├── subset_pruned_1000g.popu
      │   │   │   ├── subset_pruned_1000g.prune.in
      │   │   │   ├── subset_pruned_1000g.prune.out
      │   │   │   ├── subset_pruned_1000g.psam
      │   │   │   ├── subset_pruned_1000g.pvar
      │   │   │   ├── subset_pruned_data_1000g_key.txt
      │   │   │   ├── subset_pruned_data.log
      │   │   │   ├── subset_pruned_data_original.pvar
      │   │   │   ├── subset_pruned_data.pgen
      │   │   │   ├── subset_pruned_data.psam
      │   │   │   ├── subset_pruned_data.pvar
      │   │   │   └── subset_pruned_data_temp.pvar
      │   │   ├── freebayes_hg19.vcf
      │   │   ├── freebayes_hg19.vcf.unmap
      │   │   ├── freebayes.log
      │   │   ├── freebayes.pgen
      │   │   ├── freebayes.psam
      │   │   ├── freebayes.pvar
      │   │   ├── freebayes.pvar_original
      │   │   ├── freebayes_tmp.pvar
      │   │   ├── freebayes.vcf
      │   │   ├── pca_projection
      │   │   │   ├── final_subset_pruned_data_pcs.log
      │   │   │   ├── final_subset_pruned_data_pcs.sscore
      │   │   │   ├── subset_pruned_1000g_pcs.acount
      │   │   │   ├── subset_pruned_1000g_pcs.eigenval
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec.allele
      │   │   │   ├── subset_pruned_1000g_pcs.log
      │   │   │   ├── subset_pruned_1000g_pcs_projected.log
      │   │   │   └── subset_pruned_1000g_pcs_projected.sscore
      │   │   └── pca_sex_checks_original
      │   │       ├── ancestry_assignments.tsv
      │   │       ├── Ancestry_PCAs.png
      │   │       └── variables.tsv
      │   ├── individual_2
      │   │   ├── common_snps
      │   │   │   ├── final_subset_pruned_data.bed
      │   │   │   ├── final_subset_pruned_data.bim
      │   │   │   ├── final_subset_pruned_data.fam
      │   │   │   ├── final_subset_pruned_data.log
      │   │   │   ├── final_subset_pruned_data.pgen
      │   │   │   ├── final_subset_pruned_data.psam
      │   │   │   ├── final_subset_pruned_data.pvar
      │   │   │   ├── snps_1000g.tsv
      │   │   │   ├── SNPs2keep.txt
      │   │   │   ├── snps_data.tsv
      │   │   │   ├── subset_1000g.log
      │   │   │   ├── subset_1000g.pgen
      │   │   │   ├── subset_1000g.psam
      │   │   │   ├── subset_1000g.pvar
      │   │   │   ├── subset_data.log
      │   │   │   ├── subset_data.pgen
      │   │   │   ├── subset_data.prune.out
      │   │   │   ├── subset_data.psam
      │   │   │   ├── subset_data.pvar
      │   │   │   ├── subset_pruned_1000g.bed
      │   │   │   ├── subset_pruned_1000g.bim
      │   │   │   ├── subset_pruned_1000g.fam
      │   │   │   ├── subset_pruned_1000g.log
      │   │   │   ├── subset_pruned_1000g.pgen
      │   │   │   ├── subset_pruned_1000g.popu
      │   │   │   ├── subset_pruned_1000g.prune.in
      │   │   │   ├── subset_pruned_1000g.prune.out
      │   │   │   ├── subset_pruned_1000g.psam
      │   │   │   ├── subset_pruned_1000g.pvar
      │   │   │   ├── subset_pruned_data_1000g_key.txt
      │   │   │   ├── subset_pruned_data.log
      │   │   │   ├── subset_pruned_data_original.pvar
      │   │   │   ├── subset_pruned_data.pgen
      │   │   │   ├── subset_pruned_data.psam
      │   │   │   ├── subset_pruned_data.pvar
      │   │   │   └── subset_pruned_data_temp.pvar
      │   │   ├── freebayes_hg19.vcf
      │   │   ├── freebayes_hg19.vcf.unmap
      │   │   ├── freebayes.log
      │   │   ├── freebayes.pgen
      │   │   ├── freebayes.psam
      │   │   ├── freebayes.pvar
      │   │   ├── freebayes.pvar_original
      │   │   ├── freebayes_tmp.pvar
      │   │   ├── freebayes.vcf
      │   │   ├── pca_projection
      │   │   │   ├── final_subset_pruned_data_pcs.log
      │   │   │   ├── final_subset_pruned_data_pcs.sscore
      │   │   │   ├── subset_pruned_1000g_pcs.acount
      │   │   │   ├── subset_pruned_1000g_pcs.eigenval
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec.allele
      │   │   │   ├── subset_pruned_1000g_pcs.log
      │   │   │   ├── subset_pruned_1000g_pcs_projected.log
      │   │   │   └── subset_pruned_1000g_pcs_projected.sscore
      │   │   └── pca_sex_checks_original
      │   │       ├── ancestry_assignments.tsv
      │   │       ├── Ancestry_PCAs.png
      │   │       └── variables.tsv
      │   ├── individual_3
      │   │   ├── common_snps
      │   │   │   ├── final_subset_pruned_data.bed
      │   │   │   ├── final_subset_pruned_data.bim
      │   │   │   ├── final_subset_pruned_data.fam
      │   │   │   ├── final_subset_pruned_data.log
      │   │   │   ├── final_subset_pruned_data.pgen
      │   │   │   ├── final_subset_pruned_data.psam
      │   │   │   ├── final_subset_pruned_data.pvar
      │   │   │   ├── snps_1000g.tsv
      │   │   │   ├── SNPs2keep.txt
      │   │   │   ├── snps_data.tsv
      │   │   │   ├── subset_1000g.log
      │   │   │   ├── subset_1000g.pgen
      │   │   │   ├── subset_1000g.psam
      │   │   │   ├── subset_1000g.pvar
      │   │   │   ├── subset_data.log
      │   │   │   ├── subset_data.pgen
      │   │   │   ├── subset_data.prune.out
      │   │   │   ├── subset_data.psam
      │   │   │   ├── subset_data.pvar
      │   │   │   ├── subset_pruned_1000g.bed
      │   │   │   ├── subset_pruned_1000g.bim
      │   │   │   ├── subset_pruned_1000g.fam
      │   │   │   ├── subset_pruned_1000g.log
      │   │   │   ├── subset_pruned_1000g.pgen
      │   │   │   ├── subset_pruned_1000g.popu
      │   │   │   ├── subset_pruned_1000g.prune.in
      │   │   │   ├── subset_pruned_1000g.prune.out
      │   │   │   ├── subset_pruned_1000g.psam
      │   │   │   ├── subset_pruned_1000g.pvar
      │   │   │   ├── subset_pruned_data_1000g_key.txt
      │   │   │   ├── subset_pruned_data.log
      │   │   │   ├── subset_pruned_data_original.pvar
      │   │   │   ├── subset_pruned_data.pgen
      │   │   │   ├── subset_pruned_data.psam
      │   │   │   ├── subset_pruned_data.pvar
      │   │   │   └── subset_pruned_data_temp.pvar
      │   │   ├── freebayes_hg19.vcf
      │   │   ├── freebayes_hg19.vcf.unmap
      │   │   ├── freebayes.log
      │   │   ├── freebayes.pgen
      │   │   ├── freebayes.psam
      │   │   ├── freebayes.pvar
      │   │   ├── freebayes.pvar_original
      │   │   ├── freebayes_tmp.pvar
      │   │   ├── freebayes.vcf
      │   │   ├── pca_projection
      │   │   │   ├── final_subset_pruned_data_pcs.log
      │   │   │   ├── final_subset_pruned_data_pcs.sscore
      │   │   │   ├── subset_pruned_1000g_pcs.acount
      │   │   │   ├── subset_pruned_1000g_pcs.eigenval
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec.allele
      │   │   │   ├── subset_pruned_1000g_pcs.log
      │   │   │   ├── subset_pruned_1000g_pcs_projected.log
      │   │   │   └── subset_pruned_1000g_pcs_projected.sscore
      │   │   └── pca_sex_checks_original
      │   │       ├── ancestry_assignments.tsv
      │   │       ├── Ancestry_PCAs.png
      │   │       └── variables.tsv
      │   ├── individual_4
      │   │   ├── common_snps
      │   │   │   ├── final_subset_pruned_data.bed
      │   │   │   ├── final_subset_pruned_data.bim
      │   │   │   ├── final_subset_pruned_data.fam
      │   │   │   ├── final_subset_pruned_data.log
      │   │   │   ├── final_subset_pruned_data.pgen
      │   │   │   ├── final_subset_pruned_data.psam
      │   │   │   ├── final_subset_pruned_data.pvar
      │   │   │   ├── snps_1000g.tsv
      │   │   │   ├── SNPs2keep.txt
      │   │   │   ├── snps_data.tsv
      │   │   │   ├── subset_1000g.log
      │   │   │   ├── subset_1000g.pgen
      │   │   │   ├── subset_1000g.psam
      │   │   │   ├── subset_1000g.pvar
      │   │   │   ├── subset_data.log
      │   │   │   ├── subset_data.pgen
      │   │   │   ├── subset_data.prune.out
      │   │   │   ├── subset_data.psam
      │   │   │   ├── subset_data.pvar
      │   │   │   ├── subset_pruned_1000g.bed
      │   │   │   ├── subset_pruned_1000g.bim
      │   │   │   ├── subset_pruned_1000g.fam
      │   │   │   ├── subset_pruned_1000g.log
      │   │   │   ├── subset_pruned_1000g.pgen
      │   │   │   ├── subset_pruned_1000g.popu
      │   │   │   ├── subset_pruned_1000g.prune.in
      │   │   │   ├── subset_pruned_1000g.prune.out
      │   │   │   ├── subset_pruned_1000g.psam
      │   │   │   ├── subset_pruned_1000g.pvar
      │   │   │   ├── subset_pruned_data_1000g_key.txt
      │   │   │   ├── subset_pruned_data.log
      │   │   │   ├── subset_pruned_data_original.pvar
      │   │   │   ├── subset_pruned_data.pgen
      │   │   │   ├── subset_pruned_data.psam
      │   │   │   ├── subset_pruned_data.pvar
      │   │   │   └── subset_pruned_data_temp.pvar
      │   │   ├── freebayes_hg19.vcf
      │   │   ├── freebayes_hg19.vcf.unmap
      │   │   ├── freebayes.log
      │   │   ├── freebayes.pgen
      │   │   ├── freebayes.psam
      │   │   ├── freebayes.pvar
      │   │   ├── freebayes.pvar_original
      │   │   ├── freebayes_tmp.pvar
      │   │   ├── freebayes.vcf
      │   │   ├── pca_projection
      │   │   │   ├── final_subset_pruned_data_pcs.log
      │   │   │   ├── final_subset_pruned_data_pcs.sscore
      │   │   │   ├── subset_pruned_1000g_pcs.acount
      │   │   │   ├── subset_pruned_1000g_pcs.eigenval
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec.allele
      │   │   │   ├── subset_pruned_1000g_pcs.log
      │   │   │   ├── subset_pruned_1000g_pcs_projected.log
      │   │   │   └── subset_pruned_1000g_pcs_projected.sscore
      │   │   └── pca_sex_checks_original
      │   │       ├── ancestry_assignments.tsv
      │   │       ├── Ancestry_PCAs.png
      │   │       └── variables.tsv
      │   ├── individual_5
      │   │   ├── common_snps
      │   │   │   ├── final_subset_pruned_data.bed
      │   │   │   ├── final_subset_pruned_data.bim
      │   │   │   ├── final_subset_pruned_data.fam
      │   │   │   ├── final_subset_pruned_data.log
      │   │   │   ├── final_subset_pruned_data.pgen
      │   │   │   ├── final_subset_pruned_data.psam
      │   │   │   ├── final_subset_pruned_data.pvar
      │   │   │   ├── snps_1000g.tsv
      │   │   │   ├── SNPs2keep.txt
      │   │   │   ├── snps_data.tsv
      │   │   │   ├── subset_1000g.log
      │   │   │   ├── subset_1000g.pgen
      │   │   │   ├── subset_1000g.psam
      │   │   │   ├── subset_1000g.pvar
      │   │   │   ├── subset_data.log
      │   │   │   ├── subset_data.pgen
      │   │   │   ├── subset_data.prune.out
      │   │   │   ├── subset_data.psam
      │   │   │   ├── subset_data.pvar
      │   │   │   ├── subset_pruned_1000g.bed
      │   │   │   ├── subset_pruned_1000g.bim
      │   │   │   ├── subset_pruned_1000g.fam
      │   │   │   ├── subset_pruned_1000g.log
      │   │   │   ├── subset_pruned_1000g.pgen
      │   │   │   ├── subset_pruned_1000g.popu
      │   │   │   ├── subset_pruned_1000g.prune.in
      │   │   │   ├── subset_pruned_1000g.prune.out
      │   │   │   ├── subset_pruned_1000g.psam
      │   │   │   ├── subset_pruned_1000g.pvar
      │   │   │   ├── subset_pruned_data_1000g_key.txt
      │   │   │   ├── subset_pruned_data.log
      │   │   │   ├── subset_pruned_data_original.pvar
      │   │   │   ├── subset_pruned_data.pgen
      │   │   │   ├── subset_pruned_data.psam
      │   │   │   ├── subset_pruned_data.pvar
      │   │   │   └── subset_pruned_data_temp.pvar
      │   │   ├── freebayes_hg19.vcf
      │   │   ├── freebayes_hg19.vcf.unmap
      │   │   ├── freebayes.log
      │   │   ├── freebayes.pgen
      │   │   ├── freebayes.psam
      │   │   ├── freebayes.pvar
      │   │   ├── freebayes.pvar_original
      │   │   ├── freebayes_tmp.pvar
      │   │   ├── freebayes.vcf
      │   │   ├── pca_projection
      │   │   │   ├── final_subset_pruned_data_pcs.log
      │   │   │   ├── final_subset_pruned_data_pcs.sscore
      │   │   │   ├── subset_pruned_1000g_pcs.acount
      │   │   │   ├── subset_pruned_1000g_pcs.eigenval
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec
      │   │   │   ├── subset_pruned_1000g_pcs.eigenvec.allele
      │   │   │   ├── subset_pruned_1000g_pcs.log
      │   │   │   ├── subset_pruned_1000g_pcs_projected.log
      │   │   │   └── subset_pruned_1000g_pcs_projected.sscore
      │   │   └── pca_sex_checks_original
      │   │       ├── ancestry_assignments.tsv
      │   │       ├── Ancestry_PCAs.png
      │   │       └── variables.tsv
      │   └── individual_6
      │       ├── common_snps
      │       │   ├── final_subset_pruned_data.bed
      │       │   ├── final_subset_pruned_data.bim
      │       │   ├── final_subset_pruned_data.fam
      │       │   ├── final_subset_pruned_data.log
      │       │   ├── final_subset_pruned_data.pgen
      │       │   ├── final_subset_pruned_data.psam
      │       │   ├── final_subset_pruned_data.pvar
      │       │   ├── snps_1000g.tsv
      │       │   ├── SNPs2keep.txt
      │       │   ├── snps_data.tsv
      │       │   ├── subset_1000g.log
      │       │   ├── subset_1000g.pgen
      │       │   ├── subset_1000g.psam
      │       │   ├── subset_1000g.pvar
      │       │   ├── subset_data.log
      │       │   ├── subset_data.pgen
      │       │   ├── subset_data.prune.out
      │       │   ├── subset_data.psam
      │       │   ├── subset_data.pvar
      │       │   ├── subset_pruned_1000g.bed
      │       │   ├── subset_pruned_1000g.bim
      │       │   ├── subset_pruned_1000g.fam
      │       │   ├── subset_pruned_1000g.log
      │       │   ├── subset_pruned_1000g.pgen
      │       │   ├── subset_pruned_1000g.popu
      │       │   ├── subset_pruned_1000g.prune.in
      │       │   ├── subset_pruned_1000g.prune.out
      │       │   ├── subset_pruned_1000g.psam
      │       │   ├── subset_pruned_1000g.pvar
      │       │   ├── subset_pruned_data_1000g_key.txt
      │       │   ├── subset_pruned_data.log
      │       │   ├── subset_pruned_data_original.pvar
      │       │   ├── subset_pruned_data.pgen
      │       │   ├── subset_pruned_data.psam
      │       │   ├── subset_pruned_data.pvar
      │       │   └── subset_pruned_data_temp.pvar
      │       ├── freebayes_hg19.vcf
      │       ├── freebayes_hg19.vcf.unmap
      │       ├── freebayes.log
      │       ├── freebayes.pgen
      │       ├── freebayes.psam
      │       ├── freebayes.pvar
      │       ├── freebayes.pvar_original
      │       ├── freebayes_tmp.pvar
      │       ├── freebayes.vcf
      │       ├── pca_projection
      │       │   ├── final_subset_pruned_data_pcs.log
      │       │   ├── final_subset_pruned_data_pcs.sscore
      │       │   ├── subset_pruned_1000g_pcs.acount
      │       │   ├── subset_pruned_1000g_pcs.eigenval
      │       │   ├── subset_pruned_1000g_pcs.eigenvec
      │       │   ├── subset_pruned_1000g_pcs.eigenvec.allele
      │       │   ├── subset_pruned_1000g_pcs.log
      │       │   ├── subset_pruned_1000g_pcs_projected.log
      │       │   └── subset_pruned_1000g_pcs_projected.sscore
      │       └── pca_sex_checks_original
      │           ├── ancestry_assignments.tsv
      │           ├── Ancestry_PCAs.png
      │           └── variables.tsv
      └── snps_1000g_common_across_sites.tsv


Inforamtive Results
--------------------

``ancestry_assignments.tsv``: This file contains the predicted ancestries for all the samples (this will be in your output directory ``/path/to/parent/out/dir/ancestry_assignments.tsv``):

+---------+---------+-----------------+---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+------------+------------+-----------------+---------+---------+---------+---------+---------+---------------------+-----------------------------+------------------+
| FID     | IID     | PC1             | PC2           | PC3             | PC4             | PC5             | PC6             | PC7             | PC8             | PC9        | PC10       | SuperPop        | AFR     | AMR     | EAS     | EUR     | SAS     | combined_assignment | Plot                        | Final_Assignment |
+=========+=========+=================+===============+=================+=================+=================+=================+=================+=================+============+============+=================+=========+=========+=========+=========+=========+=====================+=============================+==================+
| 0       | 1       | 0.137299        | -0.107681     | -0.024958       | -0.0424577      | 0.0315921       | -0.0417027      | 0.00116216      | -0.0208361      | -0.0857521 | -0.0192496 |                 | 0       | 0       | 1       | 0       | 0       | EAS                 | Projected Data Assignments  | EAS              |
+---------+---------+-----------------+---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+------------+------------+-----------------+---------+---------+---------+---------+---------+---------------------+-----------------------------+------------------+

``Ancestry_PCAs.png``: a separate figure generated for each individual in each pool. For example, a PCA plot for individual 1 in Pool 1: ``/path/to/parent/out/dir/Pool1/individual_1/pca_sex_checks_original/Ancestry_PCAs.png``.

- This figure shows the 1000G individual locations in PC space compared to the individual. For example:

  .. figure:: ../_figures/Ancestry_PCAs.png
    :figwidth: 700px


    
- There will be an ``ancestry_assignments.tsv`` file generated for each individual in each pool and one that has all the individuals joined together in the base output directory.

  - This file has the annotations and probabilities for each pool. For example:

  +---------+---------+-----------------+---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+------------+------------+-----------------+---------+---------+---------+---------+---------+---------------------+-----------------------------+------------------+
  | FID     | IID     | PC1             | PC2           | PC3             | PC4             | PC5             | PC6             | PC7             | PC8             | PC9        | PC10       | SuperPop        | AFR     | AMR     | EAS     | EUR     | SAS     | combined_assignment | Plot                        | Final_Assignment |
  +=========+=========+=================+===============+=================+=================+=================+=================+=================+=================+============+============+=================+=========+=========+=========+=========+=========+=====================+=============================+==================+
  | 0       | 1       | 0.137299        | -0.107681     | -0.024958       | -0.0424577      | 0.0315921       | -0.0417027      | 0.00116216      | -0.0208361      | -0.0857521 | -0.0192496 |                 | 0       | 0       | 1       | 0       | 0       | EAS                 | Projected Data Assignments  | EAS              |
  +---------+---------+-----------------+---------------+-----------------+-----------------+-----------------+-----------------+-----------------+-----------------+------------+------------+-----------------+---------+---------+---------+---------+---------+---------------------+-----------------------------+------------------+

Additional Results with Reference-based Exectution
----------------------------------------------------

If you have reference SNP genotypes (*i.e.* microarray or whole exome or genome sequencing-called SNPs) and decided to estimate SNP-based ancestry, you will have additional results that compare the reference and single-cell ancestry predictions.
These will be located in:

- ``/path/to/parent/out/dir/reference``: This will contain the ancestry predictions for the reference-based ancestry predictions per individual.

- ``/path/to/parent/out/dir/ref_sc_ancestry_prediction_comparison``: This will provide results on reference vs single cell based ancestry predictions. These are the main comparison files with the two most informative ones highlighted and additional information below:

.. code-block:: bash
  :emphasize-lines: 6,7

  snp_ancestry_predictions.tsv
  reference_ancestry_numbers.png
  predicted_ancestry_numbers_correct.png
  predicted_ancestry_numbers_correct_identified.png
  statistics_heatmap.png
  predicted_numbers_statistics_heatmap_combined.png
  assignments_probabilities_w_ref.png


The ``predicted_numbers_statistics_heatmap_combined.png`` figure shows the number of individuals classified to each ancestry with the single cell data and colored by if they correctly or incorrectly match the reference-annotated SNP genotype data
and the heatmap below shows some statistical metrics for quantifying the single cell derived ancestry predictions:

.. figure:: ../_figures/predicted_numbers_statistics_heatmap_combined.png
  :figwidth: 300px


The ``assignments_probabilities_w_ref_identified.png`` figure show the probability of each sample to be classified to each of the different ancestries. 
This includes the reference-based predictions (microarray or whole exome or genome sequencing data) compared to the single cell based predictions

.. figure:: ../_figures/assignments_probabilities_w_ref_identified.png
  :figwidth: 600px