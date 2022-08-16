#!/usr/bin/env python
shell.executable('bash')

### This will only be run if there are reference SNP genotypes for individuals in the dataset

rule souporcell_correlate_indivs:
    input:
        output_dict["outdir"] + "/{pool}/souporcell/cluster_genotypes.vcf"
    output:
        output_dict["outdir"] + "/{pool}/souporcell/Genotype_ID_key.txt"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["souporcell_correlate_indivs_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["souporcell_correlate_indivs_memory"]
    threads: reference_ancestry_predictions_dict["souporcell_correlate_indivs_threads"]
    params:
        vcf = snp_dict["vcf"],
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        out = output_dict["outdir"] + "/{pool}/souporcell",
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Assign_Indiv_by_Geno.R -r {params.vcf} -c {input} -o {params.out}
        """

### convert vcf to plink
rule reference_vcf2plink:
    input:
        snp_dict["vcf"]
    output:
        bed = output_dict["outdir"] + "/reference/reference.pgen",
        bim = output_dict["outdir"] + "/reference/reference.pvar",
        fam = output_dict["outdir"] + "/reference/reference.psam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_vcf2plink_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_vcf2plink_memory"]
    threads: reference_ancestry_predictions_dict["reference_vcf2plink_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        souporcell = output_dict["outdir"] + "/reference/reference"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --vcf {input} --make-pgen --out {params.souporcell}
        """



### Pull just common SNPs between two groups ###
rule reference_common_snps:
    input:
        bed = output_dict["outdir"] + "/reference/reference.pgen",
        bim = output_dict["outdir"] + "/reference/reference.pvar",
        fam = output_dict["outdir"] + "/reference/reference.psam"
    output:
        snps_data = output_dict["outdir"] +  "/reference/common_snps/snps_data.tsv",
        snps_1000g = output_dict["outdir"] +  "/reference/common_snps/snps_1000g.tsv",
        bed = output_dict["outdir"] +  "/reference/common_snps/subset_data.pgen",
        bim = output_dict["outdir"] +  "/reference/common_snps/subset_data.pvar",
        fam = output_dict["outdir"] +  "/reference/common_snps/subset_data.psam",
        bed_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_1000g.pgen",
        bim_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_1000g.pvar",
        fam_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_1000g.psam",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_common_snps_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_common_snps_memory"]
    threads: reference_ancestry_predictions_dict["reference_common_snps_threads"]
    params:
        bim_1000 = "/opt/1000G/all_phase3_filtered.pvar",
        infile =  output_dict["outdir"] + "/reference/reference",
        infile_1000g = "/opt/1000G/all_phase3_filtered",
        out = output_dict["outdir"] +  "/reference/common_snps/subset_data",
        out_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_1000g",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.bim} {params.bim_1000} > {output.snps_1000g}
        singularity exec --bind {params.bind} {params.sif} awk 'NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {params.bim_1000} {input.bim} > {output.snps_data}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --extract {output.snps_data} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} --extract {output.snps_1000g} --make-pgen --out {params.out_1000g}
        """

### Prune with --indep,
rule reference_prune_1000g:
    input:
        bed_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_1000g.pgen",
        bim_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_1000g.pvar",
        fam_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_1000g.psam",
        bim = output_dict["outdir"] +  "/reference/common_snps/subset_data.pvar",
        bed = output_dict["outdir"] +  "/reference/common_snps/subset_data.pgen",
        fam = output_dict["outdir"] +  "/reference/common_snps/subset_data.psam",
    output:
        prune_out_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g.prune.out",
        prune_out = output_dict["outdir"] +  "/reference/common_snps/subset_data.prune.out",
        pgen_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g.pgen",
        pvar_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g.pvar",
        psam_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g.psam",
        bed_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g.bed",
        bim_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g.bim",
        fam_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g.fam",
        bed = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data.pgen",
        bim = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data.pvar",
        bim_temp = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data_temp.pvar",
        bim_old = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data_original.pvar",
        fam = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data.psam",
        data_1000g_key = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data_1000g_key.txt",
        SNPs2keep = output_dict["outdir"] +  "/reference/common_snps/SNPs2keep.txt",
        onek_popu = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g.popu"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_prune_1000g_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_prune_1000g_memory"]
    threads: reference_ancestry_predictions_dict["reference_prune_1000g_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        out_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g",
        infile_1000g = output_dict["outdir"] +  "/reference/common_snps/subset_1000g",
        infile = output_dict["outdir"] +  "/reference/common_snps/subset_data",
        out = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} \
            --indep-pairwise 50 5 0.5 \
            --out {params.out_1000g}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} --extract {output.prune_out_1000g} --make-pgen --out {params.out_1000g}
        if [[ $(grep "##" {input.bim} | wc -l) > 0 ]]
        then
            singularity exec --bind {params.bind} {params.sif} grep "##" {input.bim} > {output.data_1000g_key}
        fi
        singularity exec --bind {params.bind} {params.sif} awk -F"\\t" 'BEGIN{{OFS=FS = "\\t"}} NR==FNR{{a[$1 FS $2 FS $4 FS $5] = $0; next}} {{ind = $1 FS $2 FS $4 FS $5}} ind in a {{print a[ind], $3}}' {output.pvar_1000g} {input.bim} | singularity exec --bind {params.bind} {params.sif} grep -v "##" >> {output.data_1000g_key}
        singularity exec --bind {params.bind} {params.sif} grep -v "##" {output.data_1000g_key} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print $NF}}' > {output.prune_out}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --extract {output.prune_out} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        singularity exec --bind {params.bind} {params.sif} cp {output.bim} {output.bim_old}
        singularity exec --bind {params.bind} {params.sif} grep -v "#" {output.bim_old} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($3)}}' > {output.SNPs2keep}
        singularity exec --bind {params.bind} {params.sif} grep "#CHROM" {output.data_1000g_key} > {output.bim}
        singularity exec --bind {params.bind} {params.sif} grep -Ff {output.SNPs2keep} {output.data_1000g_key} >> {output.bim}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}NF{{NF-=1}};1' < {output.bim} > {output.bim_temp}
        singularity exec --bind {params.bind} {params.sif} grep "##" {output.pvar_1000g} > {output.bim}
        singularity exec --bind {params.bind} {params.sif} cat {output.bim_temp} >> {output.bim}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.out_1000g} --make-bed --out {params.out_1000g}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print(0,$1,$5)}}' {params.out_1000g}.psam | singularity exec --bind {params.bind} {params.sif} tail -n +2 > {params.out_1000g}.popu
        """
        
rule reference_final_pruning: ### put in contingency for duplicated snps - remove from both 1000G and your dataset
    input:
        bed = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data.pgen",
        bim = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data.pvar",
        fam = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data.psam",
    output:
        pgen = output_dict["outdir"] + "/reference/common_snps/final_subset_pruned_data.pgen",
        pvar = output_dict["outdir"] + "/reference/common_snps/final_subset_pruned_data.pvar",
        psam = output_dict["outdir"] + "/reference/common_snps/final_subset_pruned_data.psam",
        bed = output_dict["outdir"] +  "/reference/common_snps/final_subset_pruned_data.bed",
        bim = output_dict["outdir"] +  "/reference/common_snps/final_subset_pruned_data.bim",
        fam = output_dict["outdir"] +  "/reference/common_snps/final_subset_pruned_data.fam",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_final_pruning_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_final_pruning_memory"]
    threads: reference_ancestry_predictions_dict["reference_final_pruning_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        infile = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_data",
        out = output_dict["outdir"] +  "/reference/common_snps/final_subset_pruned_data"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --rm-dup 'force-first' -threads {threads} --pfile {params.infile} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 --rm-dup 'force-first' -threads {threads} --pfile {params.infile} --make-bed --out {params.out}
        """


rule reference_split_data:
    input:
        bed = output_dict["outdir"] +  "/reference/common_snps/final_subset_pruned_data.bed" 
    output:
        bed = output_dict["outdir"] +  "/reference/common_snps/split/final_subset_pruned_data.bed",
        bed_1000 = temp(output_dict["outdir"] +  "/reference/common_snps/split/subset_pruned_1000g.bed"),
        bim_1000 = temp(output_dict["outdir"] +  "/reference/common_snps/split/subset_pruned_1000g.bim"),
        fam_1000 = temp(output_dict["outdir"] +  "/reference/common_snps/split/subset_pruned_1000g.fam")
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_split_data_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_split_data_memory"]
    threads: reference_ancestry_predictions_dict["reference_split_data_threads"]
    params:
        infile = output_dict["outdir"] +  "/reference/common_snps/final_subset_pruned_data",
        outdir = output_dict["outdir"] +  "/reference/common_snps/split/",
        ref = output_dict["outdir"] +  "/reference/common_snps/subset_pruned_1000g",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} cp {params.infile}.bed {params.outdir}
        singularity exec --bind {params.bind} {params.sif} cp {params.infile}.bim {params.outdir}
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^X/23/g' {params.outdir}/final_subset_pruned_data.bim
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^Y/24/g' {params.outdir}/final_subset_pruned_data.bim
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^XY/25/g' {params.outdir}/final_subset_pruned_data.bim
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^MT/26/g' {params.outdir}/final_subset_pruned_data.bim
        singularity exec --bind {params.bind} {params.sif} cp {params.infile}.fam {params.outdir}
        singularity exec --bind {params.bind} {params.sif} cp {params.ref}.bed {params.outdir}
        singularity exec --bind {params.bind} {params.sif} cp {params.ref}.bim {params.outdir}
        singularity exec --bind {params.bind} {params.sif} cp {params.ref}.fam {params.outdir}
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^X/23/g' {params.outdir}/subset_pruned_1000g.bim
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^Y/24/g' {params.outdir}/subset_pruned_1000g.bim
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^XY/25/g' {params.outdir}/subset_pruned_1000g.bim
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^MT/26/g' {params.outdir}/subset_pruned_1000g.bim
        """


### use PCA from plink for PCA and projection
rule reference_pca_1000g:
    input:
        bed_1000g = output_dict["outdir"] + "/reference/common_snps/subset_pruned_1000g.pgen",
        bim_1000g = output_dict["outdir"] + "/reference/common_snps/subset_pruned_1000g.pvar",
        fam_1000g = output_dict["outdir"] + "/reference/common_snps/subset_pruned_1000g.psam",
        bed = output_dict["outdir"] + "/reference/common_snps/subset_pruned_data.pgen" 
    output:
        out = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs.acount",
        eig_all = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele",
        eig_vec = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs.eigenvec",
        eig = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs.eigenval",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_pca_1000g_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_pca_1000g_memory"]
    threads: reference_ancestry_predictions_dict["reference_pca_1000g_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        infile = output_dict["outdir"] + "/reference/common_snps/subset_pruned_1000g",
        out = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} \
            --freq counts \
            --pca allele-wts \
            --out {params.out}
        """


### use plink pca results to plot with R ###
rule reference_pca_project:
    input:
        bed = output_dict["outdir"] + "/reference/common_snps/final_subset_pruned_data.pgen",
        bim = output_dict["outdir"] + "/reference/common_snps/final_subset_pruned_data.pvar",
        fam = output_dict["outdir"] + "/reference/common_snps/final_subset_pruned_data.psam",
        frq = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs.acount",
        scores = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele"
    output:
        projected_scores = output_dict["outdir"] + "/reference/pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs_projected.sscore"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_pca_project_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_pca_project_memory"]
    threads: reference_ancestry_predictions_dict["reference_pca_project_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        infile = output_dict["outdir"] + "/reference/common_snps/final_subset_pruned_data",
        infile_1000g = output_dict["outdir"] + "/reference/common_snps/subset_pruned_1000g",
        out = output_dict["outdir"] + "/reference/pca_projection/final_subset_pruned_data_pcs",
        out_1000g = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs_projected"
    shell:
        """
        export OMP_NUM_THREADS={threads}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation \
                    variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} \
            --read-freq {input.frq} \
            --score {input.scores} 2 5 header-read no-mean-imputation \
                    variance-standardize \
            --score-col-nums 6-15 \
            --out {params.out_1000g}
       """


rule reference_pca_projection_assign_original:
    input:
        projected_scores = output_dict["outdir"] + "/reference/pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = output_dict["outdir"] + "/reference/pca_projection/subset_pruned_1000g_pcs_projected.sscore",
        fam_1000g = output_dict["outdir"] + "/reference/common_snps/subset_1000g.psam",
    output:
        anc_fig = output_dict["outdir"] + "/reference/pca_sex_checks_original/Ancestry_PCAs.png",
        tsv = output_dict["outdir"] + "/reference/pca_sex_checks_original/ancestry_assignments.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_pca_projection_assign_original_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_pca_projection_assign_original_memory"]
    threads: reference_ancestry_predictions_dict["reference_pca_projection_assign_original_threads"]
    params:
        variables = output_dict["outdir"] + "/reference/pca_sex_checks_original/variables.tsv",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        outdir = output_dict["outdir"] + "/reference/pca_sex_checks_original/",
        script = "/opt/ancestry_prediction_scRNAseq/scripts/PCA_Projection_Plotting_original.R",
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.outdir} > {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.projected_scores} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.projected_1000g_scores} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.fam_1000g} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {params.variables}
        """


rule reference_souporcell_comparison:
    input:
        tsv = output_dict["outdir"] + "/reference/pca_sex_checks_original/ancestry_assignments.tsv",
        sc_data = expand(output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_sex_checks_original/ancestry_assignments.tsv", zip, pool=samples.Pool, individual=samples.Individual)
    output:
        anc_fig = output_dict["outdir"] + "/ref_sc_ancestry_prediction_comparison/assignments_probabilities_w_ref_identified.png"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_souporcell_comparison_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * reference_ancestry_predictions_dict["reference_souporcell_comparison_memory"]
    threads: reference_ancestry_predictions_dict["reference_souporcell_comparison_threads"]
    params:
        variables = output_dict["outdir"] + "/ref_sc_ancestry_prediction_comparison/variables.tsv",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        outdir = output_dict["outdir"] + "/ref_sc_ancestry_prediction_comparison/",
        script = "/opt/ancestry_prediction_scRNAseq/scripts/compare_ref_seq_snp_ancestries.R",
        datadir = output_dict["outdir"],
        samples_file = input_dict["metadata_file"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.datadir} > {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.tsv} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {params.outdir} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {params.samples_file} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {params.variables}
        """
