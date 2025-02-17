#Differential 3'  differential 3' end processing of transcriptional loci using read alignments in BAM forma
rule filter_spurious_tpe_from_oversplitting:
    input:
        bam='data/aligned_data/{sample_name}.genome.bam',
        seq_sum='data/sequencing_summaries/{sample_name}_sequencing_summary.txt'
    output:
        bam='data/aligned_data/{sample_name,[^.]+}.filtered.bam'
    conda:
        'env_configs/d3pendr.yaml'
    shell:
        '''
        python scripts/filter_nanopore_oversplitting.py \
          -b {input.bam} \
          -s {input.seq_sum} \
          -o {output}
        '''
        
rule index_filtered_bam:
    input:
        bam = 'data/aligned_data/{sample_name}.filtered.bam'
    output:
        idx = 'data/aligned_data/{sample_name}.filtered.bam.bai'
    conda:
        'env_configs/samtools.yaml'
    shell:
        "samtools index {input.bam}"

def sample_name_subset(cond):
    sample_names = glob_wildcards(
        'data/raw_data/{sample_name}'
    ).sample_name
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names

def d3pendr_input(wildcards):
    return {
        'cntrl_bams': expand(
             'data/aligned_data/{sample_name}.filtered.bam',
             sample_name=sample_name_subset(wildcards.cntrl)
        ),
        'treat_bams': expand(
             'data/aligned_data/{sample_name}.filtered.bam',
             sample_name=sample_name_subset(wildcards.treat)
        ),
                'cntrl_bams.idx': expand(
             'data/aligned_data/{sample_name}.filtered.bam.bai',
             sample_name=sample_name_subset(wildcards.cntrl)
        ),
                'treat_bams.idx': expand(
             'data/aligned_data/{sample_name}.filtered.bam.bai',
             sample_name=sample_name_subset(wildcards.treat)
        )
    }

rule run_d3pendr:
    input:
        unpack(d3pendr_input),
        gtf=config['gtf_fn'],
    output:
        'data/apa_results/{treat}_vs_{cntrl}.apa.bed'
    params:
        cntrl_flag=lambda wc, input: ' '.join([f'-c {fn}' for fn in input.cntrl_bams]),
        treat_flag=lambda wc, input: ' '.join([f'-t {fn}' for fn in input.treat_bams]),
        output_prefix=lambda wc: f'data/apa/{wc.treat}_vs_{wc.cntrl}',
        bootstraps=config['d3pendr_parameters'].get('nboots', 999),
        min_read_overlap=config['d3pendr_parameters'].get('min_read_overlap', 0.2),
        use_model='--use-gamma-model' if config['d3pendr_parameters'].get('use_gamma_model', True) \
                                      else '--no-model',
    threads: 24
    conda:
        'env_configs/d3pendr.yaml'
    shell:
        '''
        d3pendr \
          {params.cntrl_flag} \
          {params.treat_flag} \
          -a {input.gtf} \
          -o {params.output_prefix} \
          --write-apa-sites \
          -p {threads} \
          --bootstraps {params.bootstraps} \
          --min-read-overlap {params.min_read_overlap} \
          {params.use_model}
        '''
