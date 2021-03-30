# HiCnv version 3

HiCnv is a pipeline to call CNVs from Hi-C data and is now using Snakemake to
facilitate the usage and interface. Currently HiCnv is setup to download 
SRA files followed by alignment with HiCPro and CNV analysis. 

# Preparing the reference files 
## hg38 reference files
Rule download_hg38_files

Rule bowtie2_index_ref_genome

Rule digest_reference_genome: digest hg38 using Mboi and HindIII

snakemake --profiles profile/pbs-torque `completed documentation coming`

# Preparing the HiCPro config
Follow instructions at within the HiCPro documentation to generate an appropriate config file
Github link (https://github.com/nservant/HiC-Pro)


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
