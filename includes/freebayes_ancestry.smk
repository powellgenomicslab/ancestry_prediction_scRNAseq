#!/usr/bin/env python
shell.executable('bash')




checkpoint subset_bam:
    input:
        bam = lambda wildcards: scrnaseq_libs_df["Bam_Files"][wildcards.pool],
        cells = lambda wildcards: scrnaseq_libs_df["Annotated_Barcode_Files"][wildcards.pool]
    output:
        output_dict["outdir"] + "/{pool}/bams/subset_bam.done"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["subset_bam_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["subset_bam_memory"]
    threads: freebayes_ancestry_dict["subset_bam_threads"]
    params:
        barcode_tag = input_dict["barcode_tag"],
        out = output_dict["outdir"] + "/{pool}/bams/",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
    shell:
        """
        /directflow/SCCGGroupShare/projects/DrewNeavin/software/anaconda3/envs/sinto/bin/sinto filterbarcodes -b {input.bam} -c {input.cells} --barcodetag {params.barcode_tag} --outdir {params.out} --nproc {threads}
        echo "done" > {output}
        """
        # singularity exec --bind {params.bind} {params.sif} sinto filterbarcodes -b {input.bam} -c {input.cells} --barcodetag {params.barcode_tag} --outdir {params.out} --nproc {threads}



rule index:
    input:
        output_dict["outdir"] + "/{pool}/bams/subset_bam.done"
    output:
        output_dict["outdir"] + "/{pool}/bams/{individual}.bam.bai"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["index_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["index_memory"]
    threads: freebayes_ancestry_dict["index_threads"]
    params:
        bam = output_dict["outdir"] + "/{pool}/bams/{individual}.bam",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} samtools index {params.bam}
        """



rule parallel_freebayes_regions:
    input:
        input_dict["metadata_file"]
    output:
        regions = output_dict["outdir"] + "/freebayes_regions_file"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["parallel_freebayes_regions_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["parallel_freebayes_regions_memory"]
    threads: freebayes_ancestry_dict["parallel_freebayes_regions_threads"]
    params:
        fasta = fasta,
        bind = input_dict["bind_path"],
        sif = input_dict["singularity_image"],
        regions = freebayes_ancestry_dict["parallel_freebayes_regions_N"]
    shell:
        """
        export TMPDIR=/tmp
        singularity exec --bind {params.bind},/tmp {params.sif} fasta_generate_regions.py {params.fasta}.fai {params.regions} > {output.regions}
        """


rule freebayes:
    input:
        bai = output_dict["outdir"] + "/{pool}/bams/{individual}.bam.bai",
        bam_done = output_dict["outdir"] + "/{pool}/bams/subset_bam.done",
        regions = output_dict["outdir"] + "/freebayes_regions_file"
    output:
        output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.vcf"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_memory"]
    threads: freebayes_ancestry_dict["freebayes_threads"]
    params:
        bam = output_dict["outdir"] + "/{pool}/bams/{individual}.bam",
        sif = input_dict["singularity_image"],
        fasta = fasta,
        bind = input_dict["bind_path"],
        bed = bed
    shell:
        """
        export TMPDIR=/tmp
        singularity exec --bind {params.bind},/tmp {params.sif} freebayes-parallel {input.regions} {threads} -f {params.fasta} -iXu -C 2 -q 20 -n 3 -E 1 -m 30 --min-coverage 6 --limit-coverage 100000 --targets {params.bed} {params.bam} > {output}
        """



if ref_dict["genome"] == "hg38":
    rule freebayes_update_vcf:
        input:
            freebayes = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.vcf"
        output:
            vcf = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes_updated_ids.vcf",
            vcf19 = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes_hg19.vcf",
            ids = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/updated_ids.tsv"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_update_vcf_memory"],
            disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_update_vcf_memory"]
        threads: freebayes_ancestry_dict["freebayes_update_vcf_threads"]
        params:
            sif = input_dict["singularity_image"],
            bind = input_dict["bind_path"],
            fasta = fasta_lift,
            chain_crossmap = chain_cross
        shell:
            """
            singularity exec --bind {params.bind} {params.sif} CrossMap.py vcf {params.chain_crossmap} {output.vcf} {params.fasta} {output.vcf19}
            """

else: 
    rule freebayes_update_vcf:
        input:
            freebayes = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.vcf"
        output:
            vcf = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes_updated_ids.vcf",
            vcf19 = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes_hg19.vcf",
            ids = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/updated_ids.tsv"
        resources:
            mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_update_vcf_memory"],
            disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_update_vcf_memory"]
        threads: freebayes_ancestry_dict["freebayes_update_vcf_threads"]
        params:
            sif = input_dict["singularity_image"],
            bind = input_dict["bind_path"]
        shell:
            """
            singularity exec --bind {params.bind} {params.sif} cp {output.vcf} {output.vcf19}
            """


rule freebayes_vcf2plink:
    input:
        freebayes = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes_hg19.vcf"
    output:
        bed = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.pgen",
        bim = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.pvar",
        fam = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.psam"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_vcf2plink_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_vcf2plink_memory"]
    threads: freebayes_ancestry_dict["freebayes_vcf2plink_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        freebayes = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --vcf {input.freebayes} --make-pgen --out {params.freebayes} --max-alleles 2
        singularity exec --bind {params.bind} {params.sif} cp {params.freebayes}.pvar {params.freebayes}.pvar_original
        singularity exec --bind {params.bind} {params.sif} sed -i 's/^chr//g' {params.freebayes}.pvar
        singularity exec --bind {params.bind} {params.sif} grep "#" {params.freebayes}.pvar > {params.freebayes}_tmp.pvar
        singularity exec --bind {params.bind} {params.sif} grep -v "#" {params.freebayes}.pvar | awk 'BEGIN{{FS=OFS="\\t"}}{{print $1 FS $2 FS $1 "_" $2 "_" $4 "_" $5 FS $4 FS $5 FS $6 FS $7}}' >> {params.freebayes}_tmp.pvar
        singularity exec --bind {params.bind} {params.sif} cp {params.freebayes}_tmp.pvar {params.freebayes}.pvar
        """


### Pull just common SNPs between two groups ###
rule freebayes_common_snps:
    input:
        bed = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.pgen",
        bim = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.pvar",
        fam = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.psam"
    output:
        snps_data = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/snps_data.tsv",
        snps_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/snps_1000g.tsv",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_common_snps_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_common_snps_memory"]
    threads: freebayes_ancestry_dict["freebayes_common_snps_threads"]
    params:
        bim_1000 = "/opt/1000G/all_phase3_filtered.pvar",
        infile = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes",
        infile_1000g = "/opt/1000G/all_phase3_filtered",
        out = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data",
        out_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {input.bim} {params.bim_1000} | sed '/^$/d' > {output.snps_1000g}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}NR==FNR{{a[$1,$2,$4,$5];next}} ($1,$2,$4,$5) in a{{print $3}}' {params.bim_1000} {input.bim} | sed '/^$/d' > {output.snps_data}
        """

### Get the SNPs that are common across all pools
rule common_snps_across_pools:
    input:
        snps_data = expand(output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/snps_data.tsv", zip, pool=samples.Pool, individual=samples.Individual)
    output:
        snps_data = output_dict["outdir"] +  "/common_snps_across_pools.tsv",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["common_snps_across_pools_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["common_snps_across_pools_memory"]
    threads: freebayes_ancestry_dict["common_snps_across_poolss_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        samples_file = input_dict["metadata_file"],
        # script = "/opt/ancestry_prediction_scRNAseq/scripts/common_snps.R",
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/ancestry_prediction_from_scRNA-seq/ancestry_prediction_scRNAseq/scripts/common_snps.R",
        outdir = output_dict["outdir"]
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {params.samples_file} {params.outdir}
        """

rule subset_common_snps:
    input:
        bed = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.pgen",
        bim = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.pvar",
        fam = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes.psam",
        snps_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/snps_1000g.tsv",
    output:
        bed = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data.pgen",
        bim = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data.pvar",
        fam = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data.psam",
        bed_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g.pgen",
        bim_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g.pvar",
        fam_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g.psam",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["subset_common_snps_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["subset_common_snps_memory"]
    threads: freebayes_ancestry_dict["subset_common_snps_threads"]
    params:
        snps = input_dict["common_snps"],
        snps_1000g = output_dict["outdir"] +  "/snps_1000g_common_across_sites.tsv",
        bim_1000 = "/opt/1000G/all_phase3_filtered.pvar",
        infile = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/freebayes",
        infile_1000g = "/opt/1000G/all_phase3_filtered",
        out = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data",
        out_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        # script = "/opt/ancestry_prediction_scRNAseq/scripts/subset_1000g_snps.R",
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/ancestry_prediction_from_scRNA-seq/ancestry_prediction_scRNAseq/scripts/subset_1000g_snps.R",
        outdir = output_dict["outdir"]
    shell:
        """
        ### First need to subset the 1000g snps for the new snps ###
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {params.snps} {input.snps_1000g} {params.outdir}

        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --extract {params.snps} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out} --rm-dup force-first
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile_1000g} --extract {params.snps_1000g} --make-pgen --out {params.out_1000g}
        """

### Prune with --indep,
rule freebayes_prune_1000g:
    input:
        bed_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g.pgen",
        bim_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g.pvar",
        fam_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g.psam",
        bim = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data.pvar",
        bed = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data.pgen",
        fam = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data.psam",
    output:
        prune_out_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.prune.out",
        prune_out = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data.prune.out",
        pgen_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.pgen",
        pvar_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.pvar",
        psam_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.psam",
        bed_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.bed",
        bim_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.bim",
        fam_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.fam",
        bed = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data.pgen",
        bim = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data.pvar",
        bim_temp = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data_temp.pvar",
        bim_old = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data_original.pvar",
        fam = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data.psam",
        data_1000g_key = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data_1000g_key.txt",
        SNPs2keep = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/SNPs2keep.txt",
        onek_popu = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.popu"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_prune_1000g_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_prune_1000g_memory"]
    threads: freebayes_ancestry_dict["freebayes_prune_1000g_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        out_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g",
        infile_1000g = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g",
        infile = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_data",
        out = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data"
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
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} --extract {output.prune_out} --rm-dup 'force-first' --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        singularity exec --bind {params.bind} {params.sif} cp {output.bim} {output.bim_old}
        singularity exec --bind {params.bind} {params.sif} grep -v "#" {output.bim_old} | singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print($3)}}' > {output.SNPs2keep}
        singularity exec --bind {params.bind} {params.sif} grep "#CHROM" {output.data_1000g_key} > {output.bim}
        singularity exec --bind {params.bind} {params.sif} grep -Ff {output.SNPs2keep} {output.data_1000g_key} >> {output.bim}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}NF{{NF-=1}};1' < {output.bim} > {output.bim_temp}
        singularity exec --bind {params.bind} {params.sif} grep "##" {output.pvar_1000g} > {output.bim}
        singularity exec --bind {params.bind} {params.sif} sed -i "/^$/d" {output.bim_temp}
        singularity exec --bind {params.bind} {params.sif} cat {output.bim_temp} >> {output.bim}
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.out_1000g} --make-bed --out {params.out_1000g}
        singularity exec --bind {params.bind} {params.sif} awk 'BEGIN{{FS=OFS="\t"}}{{print(0,$1,$5)}}' {params.out_1000g}.psam | singularity exec --bind {params.bind} {params.sif} tail -n +2 > {params.out_1000g}.popu
        """
        
rule freebayes_final_pruning: ### put in contingency for duplicated snps - remove from both 1000G and your dataset
    input:
        bed = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data.pgen",
        bim = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data.pvar",
        fam = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data.psam",
    output:
        pgen = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.pgen",
        pvar = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.pvar",
        psam = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.psam",
        bed = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.bed",
        bim = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.bim",
        fam = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.fam",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_final_pruning_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_final_pruning_memory"]
    threads: freebayes_ancestry_dict["freebayes_final_pruning_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        infile = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data",
        out = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --rm-dup 'force-first' --threads {threads} --pfile {params.infile} --make-pgen 'psam-cols='fid,parents,sex,phenos --out {params.out}
        singularity exec --bind {params.bind} {params.sif} plink2 --rm-dup 'force-first' --threads {threads} --pfile {params.infile} --make-bed --out {params.out}
        """


rule freebayes_split_data:
    input:
        bed = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.bed" 
    output:
        bed = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/split/final_subset_pruned_data.bed",
        bed_1000 = temp(output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/split/subset_pruned_1000g.bed"),
        bim_1000 = temp(output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/split/subset_pruned_1000g.bim"),
        fam_1000 = temp(output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/split/subset_pruned_1000g.fam")
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_split_data_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_split_data_memory"]
    threads: freebayes_ancestry_dict["freebayes_split_data_threads"]
    params:
        infile = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data",
        outdir = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/split/",
        ref = output_dict["outdir"] +  "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"]
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
rule freebayes_pca_1000g:
    input:
        bed_1000g = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.pgen",
        bim_1000g = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.pvar",
        fam_1000g = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g.psam",
        bed = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_data.pgen" 
    output:
        out = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs.acount",
        eig_all = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele",
        eig_vec = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs.eigenvec",
        eig = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs.eigenval",
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_pca_1000g_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_pca_1000g_memory"]
    threads: freebayes_ancestry_dict["freebayes_pca_1000g_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        infile = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g",
        out = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs"
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} plink2 --threads {threads} --pfile {params.infile} \
            --freq counts \
            --pca allele-wts \
            --out {params.out}
        """


### use plink pca results to plot with R ###
rule freebayes_pca_project:
    input:
        bed = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.pgen",
        bim = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.pvar",
        fam = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data.psam",
        frq = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs.acount",
        scores = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs.eigenvec.allele"
    output:
        projected_scores = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs_projected.sscore"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_pca_project_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_pca_project_memory"]
    threads: freebayes_ancestry_dict["freebayes_pca_project_threads"]
    params:
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        infile = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/final_subset_pruned_data",
        infile_1000g = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/subset_pruned_1000g",
        out = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/final_subset_pruned_data_pcs",
        out_1000g = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs_projected"
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


rule freebayes_pca_projection_assign_original:
    input:
        projected_scores = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/final_subset_pruned_data_pcs.sscore",
        projected_1000g_scores = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_projection/subset_pruned_1000g_pcs_projected.sscore",
        fam_1000g = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/common_snps/subset_1000g.psam",
    output:
        anc_fig = report(output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_sex_checks_original/Ancestry_PCAs.png", category = "PCA Cluster", caption = "ancestry_pca.rst"),
        tsv = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_sex_checks_original/ancestry_assignments.tsv"
    resources:
        mem_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_pca_projection_assign_original_memory"],
        disk_per_thread_gb=lambda wildcards, attempt: attempt * freebayes_ancestry_dict["freebayes_pca_projection_assign_original_memory"]
    threads: freebayes_ancestry_dict["freebayes_pca_projection_assign_original_threads"]
    params:
        variables = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_sex_checks_original/variables.tsv",
        sif = input_dict["singularity_image"],
        bind = input_dict["bind_path"],
        outdir = output_dict["outdir"] + "/{pool}/individual_{individual}/freebayes/pca_sex_checks_original/",
        # script = "/opt/ancestry_prediction_scRNAseq/scripts/PCA_Projection_Plotting_original.R",
        script = "/directflow/SCCGGroupShare/projects/DrewNeavin/ancestry_prediction_from_scRNA-seq/ancestry_prediction_scRNAseq/scripts/PCA_Projection_Plotting_original.R",
    shell:
        """
        singularity exec --bind {params.bind} {params.sif} echo {params.outdir} > {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.projected_scores} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.projected_1000g_scores} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} echo {input.fam_1000g} >> {params.variables}
        singularity exec --bind {params.bind} {params.sif} Rscript {params.script} {params.variables}
        """
