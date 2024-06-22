# ViralNGSBenchmarking

This repository contains supplementary data and code for a study titled "Comparative Evaluation of Bioinformatic Pipelines for Full-Length Viral Genome Assembly," published on BioRxiv (https://doi.org/10.1101/2024.03.13.584779).

## Project initialization folders
To enable easier reproducibility of the analysis, we provide the data and configuration files used for generating the three major datasets in our study (SGS-FULL, SS+NGS, and SIM). For the SS+NGS dataset, we also preliminarily publish the Sanger and next-generation sequencing results of the 49 Hungarian HIV samples used in the study. These datasets will be available in INSDC data repositories (GenBank and SRA) within one year from the point of publication in a scientific journal. Until then, the sequencing results uploaded to this repository can only be used for the replication and validation of our benchmarking analysis.

## Raw result files
* **Shankarappa.xlsx**: Diversity statistics of the Shankarappa et al. dataset
* **Vibholm.xlsx**: Diversity statistics of the Vibholm et al. dataset
* **GABM_NGS_results.csv**: Runtime and memory benchmarking results on the NGS (Ode et al.) dataset
* **GABM_SIM_results_revision.csv**: Benchmarking results using simulated HIV-1 quasispecies data
* **GABM_SIM_150.csv**: Benchmarking results using simulated HIV-1 quasispecies data with 2x150bp paired-end Illumina data
* **GABM_SGS_results.csv**: Benchmarking results using the Vibhol et al. SGS dataset
* **GABM_SSNGS_HUN_results.csv**: Benchmarking results using Hungarian Sanger and next generation sequencing results from the same samples

## Scripts
The scripts in this repository have been developed to benchmark four bioinformatic pipelines (shiver, SmaltAlign, viralngs, V-Pipe) for full-length viral genome assembly.
### Main scripts
* **gabm.sh**: Main pipeline script organizing different steps of the analysis for one sample.
* **gabm_batch.sh**: A wrapper script for gabm.sh. Initializes sample directory for gabm.sh using files in the batch directory.
* **gabm_batch.wrapper.sh**: A simple wrapper for gabm_batch.sh to run several different batch analyses after each other.

### Utilities
* **vngs.sh**: Script supplied to the viralngs docker container including all steps of genome assembly with additional runtime and memory benchmarking.
* **align_MAFFT.py**: Multiple alignment of .fasta files using MAFFT.
* **changeAmbiguousToA.py**: change IUPAC ambiguity codes in Sanger sequences due to software incompatibility before base calling using NGS data.
* **correctSangerByVariants.py**: Corrects Sanger sequences in position where NGS data strongly indicates a different base call to provide consistency between datasets.
* **multifastaDiversity.py**: Calculation of average pairwise Hamming distances in a multifasta alignment.
* **deleteGapsAlignment.py**: Deletes gap characters from multifasta alignments.
* **calculateCoordinatesByMotifs.py**: Calculates the coordinates of a genetic interval in a .fasta sequence based on a start and end motif.
* **cropGenomeByCoordinates.py**: Crops HIV-1 genomes (.fasta) based on a start and end coordinate on the HXB2 reference genome.
* **changeNBases.py**: Changes unidentified bases ('N') to a random nucleotide in a .fasta file.
* **changesToNBases.py**: Changes special IUPAC nucleotide ambiguity codes to 'N'.
* **compBamLiftover.py**: Compares .bam aligments mapped to highly similar reference sequences.
* **VCF.py**: Read in utilities for .vcf file in Python developed by Kamil Slowikowski.
* **VariantAnalysis.py**: Compares .vcf variant calls devised from highly similar reference sequences.

### Scripts for gathering results
* **gather_results.sh**: Goes through the folder structure generated by gabm_batch_wrapper.sh and calls gather_results.R in every sample directory. After that, initiates the merging of results by merge_results.R
* **gather_results.R**: Collects the benchmarking results in every sample directory.
* **merge_results.R**: Merges the tables generated by gather_results.R.

### Extra scripts
* **SGS_rename.sh**: Convenience function for renaming .fasta files.
* **ss_ngs_download.sh**: Convenience function for downloading data from NCBI Nucleotide and the Sequence Read Archive.
