"""
Author: Susheel Bhanu BUSI
Affiliation: ESB group LCSB UniLU
Date: [2021-02-17]
Run: snakemake -s Snakefile --use-conda --cores 5 -rp
Latest modification:
"""

import os
import glob
import pandas as pd

configfile:"config.yaml"
DATA_DIR=config['data_dir']
RESULTS_DIR=config['results_dir']
SAMPLES=[line.strip() for line in open("sample_list", 'r')]    # if using a sample list instead of putting them in a config file

###########
rule all:
    input:
        expand(os.path.join(RESULTS_DIR, "Stats/Cluster_{sample}.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "emboss/Cluster_{sample}_consensus.fa"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "diamond/Cluster_{sample}.tsv"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Stats/Cluster_{sample}_len_GC.txt"), sample=SAMPLES),
        expand(os.path.join(RESULTS_DIR, "Pairwise_identity/Cluster_{sample}.txt"), sample=SAMPLES)

################################
rule infoseq:
    input:
        os.path.join(DATA_DIR, "Cluster_{sample}.fa")
    output:
        os.path.join(RESULTS_DIR, "Stats/Cluster_{sample}.txt")
    conda:
        os.path.join("envs/emboss.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}.stats.log")
    message:
        "Running basic Stats, i.e. length, GC content on {wildcards.sample}"
    shell:
        "(date && infoseq {input} -outfile {output} -name -length -pgc -auto -nousa && date) &> {log}"

rule mafft:
    input:
        os.path.join(DATA_DIR, "Cluster_{sample}.fa")
    output:
        os.path.join(RESULTS_DIR, "mafft/Cluster_{sample}.msf")
    conda:
        os.path.join("envs/mafft.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}.mafft.log")
    message:
        "Running multiple sequence alignment on {wildcards.sample}"
    shell:
        "(date && mafft {input} > {output} && date) &> {log}"

rule consensus:
    input:
        rules.mafft.output
    output:
        os.path.join(RESULTS_DIR, "emboss/Cluster_{sample}_consensus.fa")
    conda:
        os.path.join("envs/emboss.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}.emboss.log")
    message:
        "Generating consensus sequence for {wildcards.sample}"
    shell:
        "(date && cons -sequence {input} -outseq {output} -name Cluster_{wildcards.sample} -auto && date) &> {log}"

rule uniprot:
    input:
        fa=rules.consensus.output,
        db=config["diamond"]["db"]
    output:
        daa=os.path.join(RESULTS_DIR, "diamond/Cluster_{sample}.daa"),
        tsv=os.path.join(RESULTS_DIR, "diamond/Cluster_{sample}.tsv")
    conda:
        os.path.join("envs/diamond.yaml")
    threads:
        config["diamond"]["threads"]
    params:
        outfmt="6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore qlen slen"
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}.uniprot.log")
    message:
        "Running UniProt trembl analyses for {wildcards.sample}"
    shell:
        "(date && "
        "daa={output.daa} && "
        "diamond blastx -q {input.fa} --db {input.db} --out {output.daa} -p {threads} --outfmt 100 && "
        "diamond view --daa ${{daa%.*}} --max-target-seqs 1 -p {threads} --outfmt {params.outfmt} --out {output.tsv} && "
        "date) &> {log}"

rule identity:
    input:
        rules.mafft.output
    output:
        os.path.join(RESULTS_DIR, "Pairwise_identity/Cluster_{sample}.txt")
    conda:
        os.path.join("envs/clustal.yaml")
    log:
        os.path.join(RESULTS_DIR, "logs/{sample}.identity.log")
    message:
        "Calculating pairwise distances from MSA for {wildcards.sample}"
    shell:
        "(date && clustalo -i {input} --distmat-out={output} --full --percent-id && date) &> {log}"

#########
# Stats #
#########
rule stats:
    input:
        rules.info.output,
        rules.identity.output
    output:
        os.path.join(RESULTS_DIR, "Stats/Cluster_{sample}_len_GC.txt")
    log:
        os.path.join(RESULTS_DIR, "logs/Cluster_{sample}.stats.log")
    message:
        "Estimating Mean, SD of sequence lengths, GC and PID (percent identity) for Cluster_{wildcards.sample}"
    run:
        # getting stats from infoseq output, i.e. length and GC
        df=pd.read_csv(input[0], header=0, delim_whitespace=True)
        stat=df.describe()
        stat=stat[['Length', '%GC']]

        # getting PID from MSA files
        pair=pd.read_csv(input[1], header=None, skiprows=1, delim_whitespace=True)
        stat['PID']=pair.describe()[1]
        
        # writing Stats to file
        stat.to_csv(output[0], sep='\t', index=True, header=True)
