# Code supporting 'Disruption of the mRNA m6A writer complex triggers autoimmunity in Arabidopsis'

## Snakemake pipelines

- nanopore_pipeline

  Pipeline used to basecall fast5 files, run QC on fastq files, align reads to transcriptome & genome, estimate poly(A) site and length, predict m6A modification and produce gene counts table.

- illumina_pipeline

  Pipeline used to QC fastq files, align reads to transcriptome & genome, produce gene counts table and calculate coverage for plots

## R notebooks
- flood_innoculation.Rmd
  
    Plot the results of two flood innoculation assays and run AVONA tests
  
- trypan_blue_staining.Rmd
  
    Plot the results of typan blue staining as quantified by ImageJ and run ANOVA tests

## Python notebooks

## Additional Scripts
