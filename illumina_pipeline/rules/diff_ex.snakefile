#Create the tx2gene table from gtf
rule tx2gene:
    input:
        gtf = config["gtf"]
    output:
        out = "data/ref/tx2gene.csv"
    params:
        dir = "data/ref"
    log:
        "logs/tx2gene.log"
    conda:
        "env_configs/R.yaml"
    shell:
        '''
        Rscript scripts/tx2geneDB.R --gtf {input.gtf} --output {output.out}
        '''

#Import quantifications and output combined table 
rule DSEq_import:
    input:
        sample_list = "data/samples.txt",
        tx2gene = "data/ref/tx2gene.csv",
        quant = expand("data/quant/{sample}/quant.sf", sample=config["samples"])
    params:
        dir = "data/quant"
    conda:
        "env_configs/R.yaml"
    output:
        gene_counts = "data/expression/gene.counts.table",
        transcript_counts = "data/expression/transcript.counts.table"
    shell:
        '''
        Rscript scripts/deseq_import.R --sample_list {input.sample_list} --tx2gene {input.tx2gene} --dir {params.dir} 
        cp {params.dir}/gene.counts.table {output.gene_counts}
        cp {params.dir}/transcript.counts.table {output.transcript_counts}
        '''
