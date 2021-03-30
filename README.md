<p align="center">
  <img src="https://raw.githubusercontent.com/joreynajr/HiCnv/master/images/vecteezy_soundwave_585767/sarmi1-03.jpg" width="400"/>
</p>

---

**HiCnv** works on contact counts at the single restriction enzyme (RE) fragment level in order to leverage Hi-C data at its highest possible and native resolution. Briefly, HiCnv first computes 1D read coverage for each RE fragment, followed by normalization for GC content, mappability and fragment length, and by smoothing using kernel density estimation (KDE). KDE smoothed counts are divided into potential CNV segments using a Hidden Markov Model (HMM), and these segments are further processed for refinement of their breakpoint coordinates (segment ends) and assignment of their CNV labels. See our [paper](https://academic.oup.com/bioinformatics/article/34/2/338/4557186) for more details.

The current version of this HiCnv uses Snakemake to facilitate deployment and improve reproducibility. Currently, HiCnv is setup to download SRA files followed by alignment with HiCPro and CNV analysis. 

# Process the hg38 reference files
To download the hg38 reference use (Rule download_hg38_files): 
<pre>
snakemake --profile workflow/profiles/local results/refs/hg38/hg38.fa.gz
</pre>

Then index those reference files using (Rule bowtie2_index_ref_genome):
<pre>
snakemake --profile workflow/profiles/local results/refs/hg38/hg38.1.bt2
</pre>


Lastly, digest the reference genome in with Mbo1 and HindIII, these will be used
later in the HiCPro configuration process (# Rule digest_reference_genome): 
<pre>
snakemake --profile workflow/profiles/local results/refs/restriction_enzymes/hg38_mboi_digestion.bed
snakemake --profile workflow/profiles/local results/refs/restriction_enzymes/hg38_hindiii_digestion.bed
</pre>

# Prepare the HiCPro config
The HiCPro configuration file must be setup with the file paths from the digested before. This 
part has to be done manually and to the specification of [HiCPro](https://github.com/nservant/HiC-Pro)
and store the configuration file within `/results/refs/hicpro/config-hicpro.{enzyme_name}.txt`.

# Download SRA paired fastq data
Uses Rule download_paired_fastq_sra:
<pre>
snakemake --profiles profile/pbs-torque `results/main/{cline}/sra/{srr}_1.fastq.gz
</pre>

# Alignment
Process your Hi-C fastq files with [HiCPro pipeline](https://github.com/nservant/HiC-Pro) (Rule hicpro_align):
<pre>
snakemake --profiles profile/pbs-torque `results/main/{cline}/{srr}/hicpro/{cline}.{srr}.ran.flag
</pre>

# Generate the coverage, GC content, mappability and fragment length information file
Download the mappability file which is used to generate the F_GC_MAP file (Rule download_hg38_mappability):
<pre>
snakemake --profiles profile/pbs-torque <<completed documentation coming>>
</pre>

Rule process_refeature: generate a restriction fragment specific file known as the
*.fragments.F_GC_MAP.bed (Fragment length, GC content and Mappability information file
file using existings commands.
<pre>
snakemake --profiles profile/pbs-torque <<completed documentation coming>>
</pre>

One dimensionalize your HiC data using pre-existing "scripts Read_coverage_generation/run_1DReadCoverage.pl" and create the *.perREfragStats file (Rule oned_read_coverage):
<pre>
snakemake --profiles profile/pbs-torque `not completed`
</pre>

# Run the CNV analysis
Finally run the CNV analyis which uses scripts/hicnv_v2.R (Rule run_hicnv):
<pre>
snakemake --profiles profile/pbs-torque `not completed`
</pre>

# Workflow schema
The following dag diagram was generated by snakemake and shows the connections between snakemake rules which for the basis for the HiCnv workflow.
<p align="center">
  <img src="https://raw.githubusercontent.com/joreynajr/HiCnv/master/images/graph.svg" width="400"/>
</p>


# Contact

mod: jreyna@lji.org (Joaquin Reyna)

orig: abhijit@lji.org (Abhijit Chakraborty)
