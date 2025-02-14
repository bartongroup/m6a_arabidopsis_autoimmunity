rule make gentrome:
    input:
        trans = "data/transcriptome/stringtie.fa",
        genome = config["genome"]
    output:
        gentrome = "data/transcriptome/gentrome.fa"
    shell:
        "cat {input.trans} {input.genome} > {output.gentrome}"

#Remove duplicates is currently sending a lot of unwanted output to files
rule remove_dup:
    input:
        gentrome = "data/transcriptome/gentrome.fa"
    output:
        clean = "data/transcriptome/gentrome.clean.fasta",
        dup = "data/transcriptome/gentrome.dup.fasta"
    conda:
        "env_configs/seqkit.yaml"
    shell:
        "seqkit rmdup -s -i {input.gentrome} -o {output.clean} --dup-seqs-file {output.dup}"
        
rule salmon_index:
    input:
        gentrome = "data/transcriptome/gentrome.clean.fasta",
        decoys = "data/ref/decoys.txt"
    output:
        dir = directory("data/quant/salmon_index")
    conda:
        "env_configs/salmon.yaml"
    shell:
        "salmon index -t {input.gentrome} -d {input.decoys} -i {output.dir} --gencode"
    
rule salmon_quant:
    input:
        R1 = "data/reads/{sample}_R1_001.fastq.gz",
        R2 = "data/reads/{sample}_R2_001.fastq.gz",
        index = "data/quant/salmon_index"
    output:
        quant = "data/quant/{sample}/quant.sf"
    conda:
        'env_configs/salmon.yaml'
    threads: 
        config["ncore"]
    params: 
        prefix=lambda wc: f"data/quant/{wc.sample}"
    shell:
        '''
        salmon quant -l A -p {threads} \
          --validateMappings \
          --numGibbsSamples 100 \
          --dumpEq \
          --posBias \
          -i {input.index} \
          -1 {input.R1} \
          -2 {input.R2} \
          -o {params.prefix}
        '''
        
rule get_table:
    input:
        sample_list = "data/samples.txt"
    params:
        dir = "data/quant"
    conda:
        "env_configs/R.yaml"
    output:
        transcript_counts = "data/quant/transcript.counts.table"
    shell:
        '''
        Rscript scripts/tximport.R --sample_list {input.sample_list} --dir {params.dir} 
        '''
