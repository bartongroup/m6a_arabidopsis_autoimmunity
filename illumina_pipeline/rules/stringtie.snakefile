#Use stringtie to build reference transcriptome

rule stringtie: 
    input:
        bam = "data/alignments/merged_{cond}.bam",
        gtf = config["gtf"]
    output:
        gtf = "data/transcriptome/{cond}.set.gtf"
    log:
        "logs/{cond}.stringtie.log"
    conda:
        'env_configs/stringtie.yaml'
    shell:
        "stringtie -o {output.gtf} -G {input.gtf} {input.bam}"
        
rule stringtie_merge:
    input:
        gtf_list = "data/transcriptome/gtf_list.txt",
        ref_gtf = config["gtf"],
        gtf = expand("data/transcriptome/{cond}.set.gtf", cond = config["conds"])
    output:
        gtf = "data/transcriptome/stringtie_merged.tmp.gtf"
    conda:
        'env_configs/stringtie.yaml'
    shell:
        "stringtie --merge -G {input.ref_gtf} -p 8 -o {output.gtf} {input.gtf_list}"
        
#Call on Matt's script to tidy up the gene id names produced by stringtie
rule tidy_gtf:
    input:
        gtf = "data/transcriptome/stringtie_merged.tmp.gtf"
    output:
        gtf = "data/transcriptome/stringtie_merged.gtf"
    shell:
        '''
        python scripts/relabel_stringtie_merge_gene_ids.py {input.gtf} {output.gtf}
        rm {input.gtf}
        '''
        
rule getfasta:
    input:
        gtf = "data/transcriptome/stringtie_merged.gtf"
    output:
        fasta = "data/transcriptome/stringtie.fa"
    conda:
        "env_configs/stringtie.yaml"
    params:
        ref = config["genome"]
    shell:
        "gffread -w {output.fasta} -g {params.ref} {input.gtf}"
