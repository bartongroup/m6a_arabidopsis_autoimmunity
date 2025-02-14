import os
import re

#Create a STAR index using ref genome and annotation specified in config file
rule star_index:
    input:
        ref_genome = config["genome"],
        gtf = config["gtf"]
    params:
        n_core = config["ncore"]
    output:
        directory('data/alignments/STAR_index')
    conda:
        'env_configs/star.yaml'
    shell:
        '''
        mkdir {output};
        STAR \
          --runThreadN {params.n_core} \
          --runMode genomeGenerate \
          --genomeDir {output} \
          --genomeFastaFiles {input.ref_genome} \
          --sjdbGTFfile {input.gtf} 
        '''    
        
def bam_output(wildcards):
    samples = glob_wildcards(
        'data/reads/{sample}_R1_001.fastq.gz'
    ).sample

    return expand(
        ['data/alignments/{sample}_Aligned.sortedByCoord.out.bam'],
        sample=samples
    )
        
rule star_align:
    input:
        R1 = "data/reads/{sample}_R1_001.fastq.gz",
        R2 = "data/reads/{sample}_R2_001.fastq.gz",
        index = "data/alignments/STAR_index"
    params:
        n_core = config["ncore"],
        prefix=lambda wc: f"data/alignments/{wc.sample}_"
    output:
        "data/alignments/{sample}_Aligned.sortedByCoord.out.bam"
    conda:
        'env_configs/star.yaml'
    shell:
        '''
        STAR \
          --runThreadN {params.n_core} \
          --genomeChrBinNbits 12 \
          --genomeDir {input.index} \
          --readFilesIn {input.R1} {input.R2} \
          --readFilesCommand zcat \
          --outSAMtype BAM SortedByCoordinate \
          --quantMode GeneCounts \
          --outReadsUnmapped Fastx\
          --outSAMstrandField intronMotif \
          --alignIntronMin 60 \
          --alignIntronMax 20000 \
          --chimOutType Junctions \
          --chimSegmentMin 15 \
          --outFileNamePrefix {params.prefix}
        '''
        
rule unmapped:
    input:
        bamfile = 'data/alignments/{sample}_Aligned.sortedByCoord.out.bam'
    output:
        unmapped = 'data/alignments/{sample}.unmapped'
    params:
        bam = lambda wc: f"data/alignments/{wc.sample}_Aligned.sortedByCoord.out.bam"
    conda:
        'env_configs/samtools.yaml'
    shell:
        "samtools view -hf 4 {params.bam} > {output.unmapped}"
       
rule stats:
    input:
        bam = "data/alignments/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        stats = "data/alignments/{sample}_Aligned.sortedByCoord.out.bamstats"
    conda:
        'env_configs/samtools.yaml'
    shell:
        '''
        samtools flagstat {input.bam} > {output.stats}
        '''
        
rule index_bam:
    input:
        bam = "data/alignments/{sample}_Aligned.sortedByCoord.out.bam"
    output:
        index = "data/alignments/{sample}_Aligned.sortedByCoord.out.bam.bai"
    conda:
        'env_configs/samtools.yaml'
    shell:
        "samtools index {input.bam}"

#Split bam files into condition groups
def sample_name_subset(cond):
    sample_names = glob_wildcards('data/reads/{sample}_R1_001.fastq.gz').sample
    return [sn for sn in sample if sn.rsplit('_', 1)[0] == cond]
    
rule multibam:
    input:
        bam = lambda wc: expand(
            "data/alignments/{sample}_Aligned.sortedByCoord.out.bam",
            sample=sample_name_subset(wc.cond)
        )
    output:
        cond_bam = "data/alignments/merged_{cond}.bam"
    conda:
        'env_configs/samtools.yaml'
    shell:
        "samtools merge {output.cond_bam} {input.bam}"

rule index_multibam:
    input:
        bam = "data/alignments/merged_{cond}.bam"
    output:
        index = "data/alignments/merged_{cond}.bam.bai"
    conda:
        'env_configs/samtools.yaml'
    shell:
        "samtools index {input.bam}"
