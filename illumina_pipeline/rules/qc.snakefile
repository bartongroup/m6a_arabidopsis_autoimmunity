import os
import re

rule fastqc:
    input:
        reads = "data/reads/{sample}.fastq.gz"
    output:
        qc = "data/fastqc/{sample}_fastqc.html"
    params:
        outputpath = "data/fastqc"
    conda:
        'env_configs/qc.yaml'
    shell:
        "fastqc -t $(({config[ncore]}+0)) -o {params.outputpath} {input.reads}"
        
def multiqc_input(wildcards):
    samples = glob_wildcards(
        'data/reads/{sample}.fastq.gz'
    ).sample

    return expand(
        ['data/fastqc/{sample}_fastqc.html'],
        sample=samples
    )

rule multiqc:
    input:
        multiqc_input
    output:
        report = "data/fastqc/report_quality_control.html"
    conda:
        'env_configs/qc.yaml'
    params:
        path = "data/fastqc"
    shell:
        "multiqc {params.path} --filename {output.report}"
