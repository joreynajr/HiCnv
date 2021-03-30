<p align="center">
  <img src="https://raw.githubusercontent.com/slundberg/HiCnv/images/vecteezy_soundwave_585767/sarmi1-03.jpg" width="800" />
</p>

# HiCnv version 3

HiCnv is a pipeline to call CNVs from Hi-C data and is now using Snakemake to
facilitate the usage and interface. Currently HiCnv is setup to download 
SRA files followed by alignment with HiCPro and CNV analysis. 

# Process the hg38 reference files
To download the hg38 reference use: 
`
# Rule download_hg38_files
snakemake --profile workflow/profiles/local results/refs/hg38/hg38.fa.gz
`
Then index those reference files using:
`
Rule bowtie2_index_ref_genome
snakemake --profile workflow/profiles/local results/refs/hg38/hg38.1.bt2
`

Lastly, digest the reference genome in with Mbo1 and HindIII, these will be used
later in the HiCPro configuration process.
`
# Rule digest_reference_genome
snakemake --profile workflow/profiles/local results/refs/restriction_enzymes/hg38_mboi_digestion.bed
snakemake --profile workflow/profiles/local results/refs/restriction_enzymes/hg38_hindiii_digestion.bed

`

# Prepare the HiCPro config
The HiCPro configuration file must be setup with the file paths from the digested before. This 
part has to be done manually and to the specification of HiCPro so please follow this link (https://github.com/nservant/HiC-Pro)
and store the configuration file within `/results/refs/hicpro/config-hicpro.<<enzyme-name>>.txt`.

# Download SRA paired fastq data
Rule download_paired_fastq_sra: Download raw HiC data

snakemake --profiles profile/pbs-torque `results/main/{cline}/sra/{srr}_1.fastq.gz`

# Alignment
Rule hicpro_align: Process your Hi-C fastq files with HiCPro pipeline (https://github.com/nservant/HiC-Pro)

snakemake --profiles profile/pbs-torque `results/main/{cline}/{srr}/hicpro/{cline}.{srr}.ran.flag`

# Generate the coverage, GC content, mappability and fragment length information file

5) Rule download_hg38_mappability: download the mappability file which is used to
generate the F_GC_MAP file

snakemake --profiles profile/pbs-torque "completed documentation coming"

6) Rule process_refeature: generate a restriction fragment specific file known as the
*.fragments.F_GC_MAP.bed (Fragment length, GC content and Mappability information file
file using existings commands.

snakemake --profiles profile/pbs-torque `completed documentation coming`

7) Rule: oned_read_coverage: One dimensionalize your HiC data using pre-existing
"scripts Read_coverage_generation/run_1DReadCoverage.pl" and create the *.perREfragStats file.

snakemake --profiles profile/pbs-torque `not completed`

## Run the CNV analysis

`
# Rule run_hicnv: Run hicnv_v2.R
snakemake --profiles profile/pbs-torque `not completed`
`

# Contact

mod: jreyna@lji.org (Joaquin Reyna)

orig: abhijit@lji.org (Abhijit Chakraborty)
