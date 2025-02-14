# Code supporting 'Disruption of the mRNA m6A writer complex triggers autoimmunity in Arabidopsis'

## Snakemake pipelines

- nanopore_pipeline

  Pipeline used to basecall fast5 files, run QC on fastq files, align reads to transcriptome & genome, estimate poly(A) site and length, predict m6A modification and produce gene counts table.

- illumina_pipeline

  Pipeline used to QC fastq files, align reads to transcriptome & genome, produce gene counts table and calculate coverage for plots

## R notebooks

- pr1_qpcr_plots.Rmd

- flood_innoculation.Rmd
  
    Plot the results of two flood innoculation assays and run AVONA tests
  
- trypan_blue_staining.Rmd
  
    Plot the results of typan blue staining as quantified by ImageJ and run ANOVA tests

- LCMSMS.Rmd

- m6a_stoichiometry.Rmd

  Plot the results of m6a modification prediction using m6anet and yanocomp

- polya_site_choice.Rmd

  Plot the shift in poly(A) start site

- polya_tail_lengths_temperature_m6a.Rmd

- polya_tail_difference.Rmd

- enolase_polyA.Rmd

  Plot the distribution of poly(A) tail lengths for ENOLASE II spike in sequences

## Python notebooks

- Gene track plots

## Additional Scripts
