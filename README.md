<p align="center">
  <img src="https://raw.githubusercontent.com/joreynajr/HiCnv/master/images/vecteezy_soundwave_585767/HiCnv_Logo.svg"/>
</p>

---

# HiCnv

### An optimized and flexible workflow for copy number calling from HiC data

[![DOI](https://img.shields.io/badge/DOI-10.1186%2Fs13059--015--0831--x-lightgrey.svg?style=flat-square)](https://pubmed.ncbi.nlm.nih.gov/29048467/)

----

**HiCnv** works on contact counts at the single restriction enzyme (RE) fragment level in order to leverage Hi-C data at its highest possible and native resolution. Briefly, HiCnv first computes 1D read coverage for each RE fragment, followed by normalization for GC content, mappability and fragment length, and by smoothing using kernel density estimation (KDE). KDE smoothed counts are divided into potential CNV segments using a Hidden Markov Model (HMM), and these segments are further processed for refinement of their breakpoint coordinates (segment ends) and assignment of their CNV labels. See our [paper](https://academic.oup.com/bioinformatics/article/34/2/338/4557186) for more details.

The current version of this HiCnv uses [Snakemake](https://snakemake.readthedocs.io/en/stable/) to facilitate deployment and improve reproducibility. Currently, HiCnv is setup to download SRA files followed by alignment with HiCPro and CNV analysis. PLEASE READ: Snakemake uses a system of templating which I will be using in this documentation, templating refers to the use of `{variable_name}` to generalize a file path and common templates I will use include:
- `{re}` - restriction enzyme
- `{srr}` - SRR ID
- `{lib}` - LIB ID from ENCODE (or equivalent)
- `{acc}` - short for accession and is a generalization for the srr and lib wildcards (not to be confused with Encode accessions, this definition is internal for HiCnv)
- `{cline}` - cell line (can also be though of as a sample name)

**UPDATE: HiCnv now supports downloading with ENCODE, instructions are given in the HiCnv Wiki [Run with ENCODE data](https://github.com/joreynajr/HiCnv/wiki/Run-with-ENCODE-data).** 

**Please visit the Wiki Page for other in depth topics [Wiki](https://github.com/joreynajr/HiCnv/wiki).**

## Authors

mod: jreyna@lji.org (Joaquin Reyna BSc)

orig: abhijit@lji.org (Abhijit Chakraborty PhD)

orig: ferhatay@lji.org (Ferhat Ay PhD)


## Software dependencies
Rules within this workflow attempt to facilitate the installation of a few key software but there may still be some user-specific installation steps. The following is a list of required software and packages where applicable:
- Python 2.7 
- Python 3.7
  - snakemake ([link](https://snakemake.readthedocs.io/en/stable/))
  - grabseqs ([link](https://github.com/louiejtaylor/grabseqs))
  - ...
- R 
  - changepoint ([link](https://cran.r-project.org/web/packages/changepoint/changepoint.pdf))
  - bioconductor ([link](https://bioconductor.org/install/))
  - ...
- HiCPro 3.0.0 - Git download the HiCPro repository using: `git clone git@github.com:nservant/HiC-Pro.git`
- fasterq-dump ([link](https://github.com/ncbi/sra-tools/wiki/HowTo:-fasterq-dump))

## Usage
### Snakemake config
The Snakemake config file is located within `config/config.yaml`, please set the following variables (using the YAML format already provided) with paths to the correponding resources or software:
<pre>
R4: &#60;path to R 4.X.X with library dependencies installed&#62;
python_hicpro: &#60;path to Python 3.7 with package dependencies for HiC-Pro installed ([check out more details here](https://github.com/nservant/HiC-Pro#using-hic-pro-through-conda))&#62;
hicpro_dir: &#60;path to the main HiCPro installation directory&#62; (don't need complete installation with `make` just `git clone git@github.com:nservant/HiC-Pro.git`)
python2: &#60;path to Python 2.7 with package dependencies installed&#62;
</pre>

### Process the hg38 reference files
To download the hg38 reference use (Rule download_hg38_files): 
<pre>
snakemake --profile workflow/profiles/pbs-torque results/refs/hg38/hg38.fa.gz
</pre>

Then index those reference files using (Rule bowtie2_index_ref_genome):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/refs/hg38/hg38.1.bt2
</pre>


Lastly, digest the reference genome with Mbo1 and HindIII, these will be used
later in the HiCPro configuration process (Rule digest_reference_genome): 
<pre>
snakemake --profile workflow/profiles/pbs-torque results/refs/restriction_enzymes/hg38_mboi_digestion.bed
snakemake --profile workflow/profiles/pbs-torque results/refs/restriction_enzymes/hg38_hindiii_digestion.bed
</pre>

### Prepare the HiCPro config
The HiCPro configuration file must be setup manually as specified by [HiCPro](https://github.com/nservant/HiC-Pro/blob/master/doc/MANUAL.md)
and includes setting `GENOME_FRAGMENT = results/refs/restriction_enzymes/hg38_{re}_digestion.bed`. Once completed store the configuration file within `results/refs/hicpro/config-hicpro.{re}.txt`. This process should be completed for every restriction enzyme used by your dataset.

### Download SRA paired fastq data
Uses Rule download_paired_fastq_sra:
<pre>
snakemake --profile workflow/profiles/pbs-torque results/main/{cline}/sra/{srr}_1.fastq.gz
</pre>

### Align the fastq's
Process your Hi-C fastq files with [HiCPro pipeline](https://github.com/nservant/HiC-Pro) (Rule hicpro_align_only):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/main/{cline}/hicpro/bowtie_results/bwt2/{srr}/{srr}_1_hg38.bwt2merged.bam
</pre>

### Generate the coverage, GC content, mappability and fragment length information file
Download the mappability file which is used to generate the F_GC_MAP file (Rule download_hg38_mappability):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.sorted.bedGraph
</pre>

Generate a restriction fragment specific file with the format *.fragments.F_GC_MAP.bed (Fragment length, GC content and Mappability information file, uses existings commands form HiCnv v2.0) (Rule process_refeature):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/refs/restriction_enzymes/hg38_{re}_digestion.extended.fragment.gc.map.sorted.bed
</pre>

One dimensionalize your HiC data using pre-existing "scripts Read_coverage_generation/run_1DReadCoverage.pl" and create the *.perREfragStats file (Rule oned_read_coverage):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/main/{cline}/hicnv/{srr}.perREfragStats
</pre>

### Run the CNV analysis
Finally run the CNV analyis which uses scripts/hicnv_v2.R (Rule run_hicnv):
<pre>
snakemake --profile results/main/{cline}/hicnv/{cline}.{srr}_hicnv/CNV_Estimation/{cline}.{srr}.cnv.bedGraph
</pre>

The `hicnv_v2.R` script comes from the original [HiCnv Github Repo](https://github.com/ay-lab/HiCnv) and this step/rule is the most important so below is a comprehensive description of its output tree structure:

<pre>
results/main/{cline}/hicnv/tech_run/{cline}_{srr}_hicnv/
├── {cline}.{srr}.copyNumber.txt
├── gc_map_frag.bed
├── coverage.bedGraph
├── coverage.gc_map_frag.bedGraph
│
├── normalized_data
│   └── {cline}.{srr}.F_GC_MAP.NormCount.0.2_0.5_150.bedGraph
│
├── Kernel_Smoothing
│   ├── {cline}.{srr}.chr1.counts.txt
│   ├── {cline}.{srr}.chr1.kde2d_x.txt
│   ├── {cline}.{srr}.chr1.kde2d_y.txt
│   ├── {cline}.{srr}.chr1.kde2d_z.txt
│   ├── {cline}.{srr}.chr1.param.txt
│   ├── {cline}.{srr}.chr2.counts.txt
│   ├── {cline}.{srr}.chr2.kde2d_x.txt
│   ├── {cline}.{srr}.chr2.kde2d_y.txt
│   ├── {cline}.{srr}.chr2.kde2d_z.txt
│   ├── {cline}.{srr}.chr2.param.txt
│   ├── {cline}.{srr}.chr3.counts.txt
│   ├── {cline}.{srr}.chr3.kde2d_x.txt
│   ├── {cline}.{srr}.chr3.kde2d_y.txt
│   ├── {cline}.{srr}.chr3.kde2d_z.txt
│   ├── {cline}.{srr}.chr3.param.txt
│   .   
│   .   
│   .   
│   ├── {cline}.{srr}.chrX.counts.txt
│   ├── {cline}.{srr}.chrX.kde2d_x.txt
│   ├── {cline}.{srr}.chrX.kde2d_y.txt
│   ├── {cline}.{srr}.chrX.kde2d_z.txt
│   └── {cline}.{srr}.chrX.param.txt
│
└── CNV_Estimation
    ├── {cline}.{srr}.chr1.cnv.bedGraph
    ├── {cline}.{srr}.chr2.cnv.bedGraph
    ├── {cline}.{srr}.chr3.cnv.bedGraph
    .   
    .   
    .   
    ├── {cline}.{srr}.chrX.cnv.bedGraph
    ├── {cline}.{srr}.cnv.bedGraph
    └── {cline}.{srr}.cnv.txt
</pre>

For each cell line {cline} and its technical replicates {srr} you will have the tree structure above

- {cline}.{srr}.copyNumber.txt - contains copyNumber per chromosome

- gc_map_frag.bed - contains GC and mappability information for each RE cutsite

- coverage.bedGraph - contains the coverage for each  RE cutsite

- coverage.gc_map_frag.bedGraph - contains the bedtools intersection between coverage.bedGraph and gc_map_frag.bed

- normalized_data - contains a single normalized data file 

- Kernel_Smoothing - after normalization CNV's are analyzed chromosome by chromosome. For a given chromosome 5 files are generated 
  - {cline}.{srr}.chr{num}.counts.txt - contains loci smooth count information
  - {cline}.{srr}.chr{num}.kde2d_x.txt
  - {cline}.{srr}.chr{num}.kde2d_y.txt
  - {cline}.{srr}.chr{num}.kde2d_z.txt
  - {cline}.{srr}.chr{num}.param.txt

- CNV_Estimation - similarly to the Kernel_Smoothing output, the CNV_Estimation output is done chromosome by chromosome and for a given chromosome a single file is generated, there are also two extra files which are summarizing the genome wide output:
  - {cline}.{srr}.chr{num}.cnv.bedGraph - per chromosome contains the final segmentation copy number calls
  - {cline}.{srr}.cnv.bedGraph - concatination of all chromosomes
  - {cline}.{srr}.cnv.txt

## Run the HiCnv workflow the _fast way_
After installing all the necessary software, and setting up the configurations files you can simply run (Rule run_hicnv):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/main/{cline}/hicnv/{cline}.{srr}_hicnv/CNV_Estimation/{cline}.{srr}.cnv.bedGraph
</pre>
Snakemake deduces what rules need to be run and will run as if you you followed the slower tutorial from above.

## Workflow schema
The following dag diagram was generated by snakemake and shows the connections between snakemake rules which for the basis for the HiCnv workflow.
<p align="center">
  <img src="https://raw.githubusercontent.com/joreynajr/HiCnv/master/images/graph.svg" width="800"/>
</p>

## References
- Identification of copy number variations and translocations in cancer cells from Hi-C data - https://doi.org/10.1093/bioinformatics/btx664
- Optimal Detection of Changepoints With a Linear Computational Cost - (PELT Algorithm) - https://doi.org/10.1080/01621459.2012.737745

