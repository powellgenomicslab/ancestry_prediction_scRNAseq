
#!/usr/bin/env python


##### Originally written by wheaton5 for souporcell #####
##### Edited by drneavin on 10 August, 2022 for ancestry calling pipeline #####



import argparse

parser = argparse.ArgumentParser(
    description="single cell RNAseq mixed genotype clustering using sparse mixture model clustering.")
parser.add_argument("-i", "--bam", required = True, help = "cellranger bam")
parser.add_argument("-f", "--fasta", required = True, help = "reference fasta file")
parser.add_argument("-t", "--threads", required = True, type = int, help = "max threads to use")
parser.add_argument("-o", "--out_dir", required = True, help = "name of directory to place souporcell files")
parser.add_argument("--min_alt", required = False, default = "10", help = "min alt to use locus, default = 10.")
parser.add_argument("--min_ref", required = False, default = "10", help = "min ref to use locus, default = 10.")
parser.add_argument("--common_variants", required = False, default = None, 
    help = "common variant loci or known variant loci vcf, must be vs same reference fasta")
args = parser.parse_args()



print("checking modules")
# importing all reqs to make sure things are installed
import numpy as np
import scipy
import gzip
import math
import pystan
import vcf
import pysam
import pyfaidx
import subprocess
import time
import os
print("imports done")


open_function = lambda f: gzip.open(f,"rt") if f[-3:] == ".gz" else open(f)




#test bam load
bam = pysam.AlignmentFile(args.bam)
num_cb = 0
num_cb_cb = 0 # num reads with barcodes from barcodes.tsv file
num_umi = 0
num_read_test = 100000
for (index,read) in enumerate(bam):
    if index >= num_read_test:
        break


print("checking fasta")
fasta = pyfaidx.Fasta(args.fasta, key_function = lambda key: key.split()[0])



def get_fasta_regions(fastaname, threads):
    fasta = pyfaidx.Fasta(args.fasta, key_function = lambda key: key.split()[0])
    total_reference_length = 0
    for chrom in sorted(fasta.keys()):
       total_reference_length += len(fasta[chrom])
    step_length = int(math.ceil(total_reference_length/threads))
    regions = []
    region = []
    region_so_far = 0
    chrom_so_far = 0
    for chrom in sorted(fasta.keys()):
        chrom_length = len(fasta[chrom])
        chrom_so_far = 0
        if chrom_length < 250000:
            continue
        while True:
            if region_so_far + (chrom_length - chrom_so_far) < step_length:
                region.append((chrom, chrom_so_far, chrom_length))
                region_so_far += chrom_length - chrom_so_far
                chrom_so_far = 0
                break
            else:
                region.append((chrom, chrom_so_far, chrom_so_far + step_length - region_so_far))
                regions.append(region)
                region = []
                chrom_so_far += step_length - region_so_far
                region_so_far = 0
    if len(region) > 0:
        if len(regions) == args.threads:
            regions[-1] = regions[-1] + region
        else:
            regions.append(region)
    return regions


def get_bam_regions(bamname, threads):
    bam = pysam.AlignmentFile(bamname)
    total_reference_length = 0
    for chrom in bam.references:
        total_reference_length += bam.get_reference_length(chrom)
    step_length = int(math.ceil(total_reference_length / threads))
    regions = []
    region = []
    region_so_far = 0
    chrom_so_far = 0
    for chrom in bam.references:
        chrom_length = bam.get_reference_length(chrom)
        #print(chrom+" size "+str(chrom_length)+" and step size "+str(step_length))
        while True:
            #print("\tregion so far "+str(region_so_far)+" chrom so far "+str(chrom_so_far)) 
            #print("\t"+str(chrom_length - chrom_so_far)+" <= "+str(step_length - region_so_far))
            #print("\t"+str((chrom_length - chrom_so_far) <= step_length - region_so_far))
            if (chrom_length - chrom_so_far) <= step_length - region_so_far:
                region.append((chrom, chrom_so_far, chrom_length))
                #print("\t\tending chrom\t"+chrom+":"+str(chrom_so_far)+"-"+str(chrom_length))
                region_so_far += chrom_length - chrom_so_far
                chrom_so_far = 0
                break
            else:
                region.append((chrom, chrom_so_far, chrom_so_far + step_length - region_so_far))
                #print("\t\tending region\t"+chrom+":"+str(chrom_so_far)+"-"+str(chrom_so_far + step_length - region_so_far))
                regions.append(region)
                region = []
                chrom_so_far += step_length - region_so_far
                region_so_far = 0
    if len(region) > 0:
        regions.append(region)
    
    return regions



def freebayes(args, bam, fasta):

        
    # parallelize the samtools depth call. It takes too long
    regions = get_bam_regions(bam, int(args.threads))
    depth_files = []
    depth_procs = []
    print(len(regions))
    for (index, region) in enumerate(regions):
        region_args = []
        for (chrom, start, stop) in region:
            region_args.append(chrom+":"+str(start)+"-"+str(stop))
        depthfile = args.out_dir+"/depth_"+str(index)+".bed"
        depth_files.append(depthfile)
        min_cov = int(args.min_ref)+int(args.min_alt)
        with open(depthfile, 'w') as bed:
            with open(depthfile+".sh",'w') as depther:
                depther.write("samtools view -hb "+bam+" "+" ".join(region_args)+ " | samtools depth - | "+
                "awk '{ if ($3 >= "+str(min_cov)+ " && $3 < 100000) { print $1 \"\t\" $2 \"\t\" $2+1 \"\t\" $3 } }'")
            subprocess.check_call(["chmod", "777", depthfile+".sh"])
            #ps0 = subprocess.Popen(['samtools', 'view', bam]+region_args, stdout = subprocess.PIPE)
            #ps1 = subprocess.Popen(['samtools', 'depth', '-'], stdin = ps0.stdout, stdout = subprocess.PIPE)
            # awk magic 
            #ps2 = subprocess.Popen(["awk '{ if ($3 >= " + str(min_cov) + " && $3 < 100000) { print $1 \"\t\" $2 \"\t\" $2+1 \"\t\" $3 } }'"], 
            #    shell = True, stdin = ps1.stdout, stdout = bed)
            ps = subprocess.Popen([depthfile+".sh"], shell = True, stdout = bed)
            depth_procs.append(ps)

    for proc in depth_procs:
        proc.wait()
    merged_depthfiles = []
    for depth_file in depth_files:
        merged_depthfile = depth_file[:-4]+"_merged.bed"
        with open(merged_depthfile, 'w') as bed:
            subprocess.check_call(["bedtools", "merge", "-i", depth_file], stdout = bed)
        merged_depthfiles.append(merged_depthfile)
    with open(args.out_dir + "/depth_merged.bed", 'w') as merged_bed:
        subprocess.check_call(['cat']+merged_depthfiles, stdout = merged_bed)
    for tmp in depth_files: # clean up tmp bed files
        subprocess.check_call(['rm', tmp, tmp+".sh"])
    for tmp in merged_depthfiles:
        subprocess.check_call(['rm', tmp])


    with open(args.out_dir + "/common_variants_covered_tmp.bed", 'w') as bed:
        subprocess.check_call(["bedtools", "intersect", "-wa", "-a", args.common_variants, "-b", args.out_dir + "/depth_merged.bed"], stdout = bed)
    with open(args.out_dir + "/common_variants_covered_tmp.bed") as bed:
        with open(args.common_variants) as common:
            with open(args.out_dir + "/common_variants_covered.bed",'w') as out:
                for line in common:
                    if line.startswith("#"):
                        out.write(line)
                    else:
                        break
                for line in bed:
                    out.write(line)
    with open(args.out_dir + "/variants.done", 'w') as done:
        done.write(args.out_dir + "/common_variants_covered.bed" + "\n")
    # return(args.out_dir + "/common_variants_covered.bed")


    regions = get_fasta_regions(args.fasta, int(args.threads))
    print(regions)

    region_vcfs = [[] for x in range(args.threads)]
    all_vcfs = []
    bed_files = []
    procs = [None for x in range(args.threads)]
    any_running = True
    filehandles = []
    errhandles = []

    # run renamer in parallel manner
    print("running freebayes")
    while any_running:
        any_running = False
        for (index, region) in enumerate(regions):
            block = False
            if procs[index]:
                block = procs[index].poll() == None
                if block:
                    any_running = True
                else:
                    assert not(procs[index].returncode), "freebayes subprocess terminated abnormally with code " + str(procs[index].returncode)
            if len(region_vcfs[index]) == len(region):
                block = True
            if not block:
                sub_index = len(region_vcfs[index])
                chrom = region[sub_index][0]
                start = region[sub_index][1]
                end = region[sub_index][2]
                vcf_name = args.out_dir + "/souporcell_" + str(index) + "_" + str(sub_index) + ".vcf"
                filehandle = open(vcf_name, 'w')
                filehandles.append(filehandle)
                errhandle = open(vcf_name + ".err", 'w')
                errhandles.append(errhandle)
                    
                cmd = ["freebayes", "-f", args.fasta, "-iXu", "-C", "2",
                    "-q", "20", "-n", "3", "-E", "1", "-m", "30", 
                    "--min-coverage", str(int(args.min_alt)+int(args.min_ref)), "--pooled-continuous", "--skip-coverage", "100000", "--targets", str(args.out_dir + "/common_variants_covered.bed")]
                
                cmd.extend(["-r", chrom + ":" + str(start) + "-" + str(end)])
                print(" ".join(cmd))
                cmd.append(bam)
                errhandle.write(" ".join(cmd) + "\n")
                p = subprocess.Popen(cmd, stdout = filehandle, stderr = errhandle)
                all_vcfs.append(vcf_name)
                procs[index] = p
                region_vcfs[index].append(vcf_name)
                any_running = True
        time.sleep(1)
    for filehandle in filehandles:
        filehandle.close()
    for errhandle in errhandles:
        errhandle.close()
    print("merging vcfs")
    subprocess.check_call(["ls "+args.out_dir+'/*.vcf | xargs -n1 -P'+str(args.threads)+' bgzip'],shell=True)
    all_vcfs = [vcf+".gz" for vcf in all_vcfs]
    subprocess.check_call(["ls "+args.out_dir+"/*.vcf.gz | xargs -n1 -P"+str(args.threads) +" bcftools index"],shell=True)
    with open(args.out_dir + "/souporcell_merged_vcf.vcf", 'w') as vcfout:
        subprocess.check_call(["bcftools", "concat", '-a'] + all_vcfs, stdout = vcfout)
    with open(args.out_dir + "/bcftools.err", 'w') as vcferr:
        with open(args.out_dir + "/souporcell_merged_sorted_vcf.vcf", 'w') as vcfout:
            subprocess.check_call(['bcftools', 'sort', args.out_dir + "/souporcell_merged_vcf.vcf"], stdout = vcfout, stderr = vcferr)
    if not args.common_variants == None:
        with open(args.out_dir + "/common.err", 'w') as err:
            with open(args.out_dir + "/vcftmp", 'w') as out:
                subprocess.check_call(['bedtools', 'intersect', '-wa', 
                    '-a', args.out_dir + "/souporcell_merged_vcf.vcf", '-b', args.common_variants], stdout = out, stderr = err)
        subprocess.check_call(['mv', args.out_dir + "/vcftmp", args.out_dir + "/souporcell_merged_sorted_vcf.vcf"])
    subprocess.check_call(['rm', args.out_dir + '/souporcell_merged_vcf.vcf'])
    subprocess.check_call(['bgzip', args.out_dir + "/souporcell_merged_sorted_vcf.vcf"])
    final_vcf = args.out_dir + "/souporcell_merged_sorted_vcf.vcf.gz"
    subprocess.check_call(['tabix', '-p', 'vcf', final_vcf])
    for vcf in all_vcfs:
        subprocess.check_call(['rm', vcf[:-3] + ".err"])       
        subprocess.check_call(['rm', vcf +".csi"])
    subprocess.check_call(['rm'] + all_vcfs)
    if len(bed_files) > 0:
        for bed in bed_files:
            subprocess.check_call(['rm', bed + ".bed"])
        subprocess.check_call(['rm'] + bed_files)
        
    with open(args.out_dir + "/variants.done", 'w') as done:
        done.write(final_vcf + "\n")
    return(final_vcf)


bam = args.bam

final_vcf = freebayes(args, bam, fasta)