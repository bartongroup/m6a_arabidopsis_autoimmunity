rule split_strand:
    input:
        bam='data/alignments/{sample}_Aligned.sortedByCoord.out.bam',
        bai='data/alignments/{sample}_Aligned.sortedByCoord.out.bam.bai'
    output:
        bam='data/alignments/{sample}.genome.{strand}.bam',
        bai='data/alignments/{sample}.genome.{strand}.bam.bai'
    params:
        samflags_1=lambda wc: '-f 128 -F 16' if wc.strand == 'fwd' else '-f 144',
        samflags_2=lambda wc: '-f 80' if wc.strand == 'fwd' else '-f 64 -F 16'
    threads: 4
    conda:
        'env_configs/samtools.yaml'
    shell:
        '''
        samtools view -@ {threads} -b {params.samflags_1} {input.bam} > {output.bam}.1.bam
        samtools index -@ {threads} {output.bam}.1.bam
        samtools view -@ {threads} -b {params.samflags_2} {input.bam} > {output.bam}.2.bam
        samtools index -@ {threads} {output.bam}.2.bam
        samtools merge -@ {threads} {output.bam} {output.bam}.1.bam {output.bam}.2.bam
        samtools index -@ {threads} {output.bam}
        rm {output.bam}.[12].bam
        rm {output.bam}.[12].bam.bai
        '''


rule genome_coverage:
    input:
        bam='data/alignments/{sample}.genome.{strand}.bam',
        bai='data/alignments/{sample}.genome.{strand}.bam.bai'
    output:
        'data/coverage_tracks/{sample}.{strand}.bw',
    params:
        genome=config['genome'],
    conda:
        'env_configs/samtools.yaml'
    shell:
        '''
        samtools depth -d0 {input.bam} | 
          awk -v OFS='\t' '{{print $1, $2-1, $2, $3}}' > {output}.tmp.bdg
        bedGraphToBigWig {output}.tmp.bdg <(cut -f-2 {params.genome}.fai) {output}
        rm {output}.tmp.bdg
        '''


def sample_name_subset(cond):
    sample_names = glob_wildcards(
        'data/reads/{sample}_R1_001.fastq.gz'
    ).sample
    cond_sample_names = [sn for sn in sample_names if sn.startswith(cond)]
    return cond_sample_names


rule normalised_genome_coverage:
    input:
        bams=lambda wc: expand(
            'data/alignments/{sample}.genome.{strand}.bam',
            sample=sample_name_subset(wc.cond),
            strand=wc.strand,
        )
    output:
        bw='data/coverage_tracks/pooled/{cond}.cpm.{strand}.bw'
    conda:
        'env_configs/samtools.yaml'
    threads: 4
    shell:
        '''
        samtools merge -f -@ {threads} {output.bw}.tmp.bam {input.bams}
        samtools index {output.bw}.tmp.bam
        bamCoverage --normalizeUsing CPM --binSize=1 -p {threads} \
          -b {output.bw}.tmp.bam \
          -o {output.bw}
        rm {output.bw}.tmp.bam*
        '''
        
rule raw_genome_coverage:
    input:
        bams=lambda wc: expand(
            'data/alignments/{sample}.genome.{strand}.bam',
            sample=sample_name_subset(wc.cond),
            strand=wc.strand,
        )
    output:
        bw='data/coverage_tracks/pooled/{cond}.raw.{strand}.bw'
    conda:
        'env_configs/samtools.yaml'
    threads: 4
    shell:
        '''
        samtools merge -f -@ {threads} {output.bw}.tmp.bam {input.bams}
        samtools index {output.bw}.tmp.bam
        bamCoverage --normalizeUsing None --binSize=1 -p {threads} \
          -b {output.bw}.tmp.bam \
          -o {output.bw}
        rm {output.bw}.tmp.bam*
        '''
