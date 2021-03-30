<p align="center">
  <img src="https://raw.githubusercontent.com/joreynajr/HiCnv/master/images/vecteezy_soundwave_585767/HiCnv_Logo.svg"/>
</p>

---

**HiCnv** works on contact counts at the single restriction enzyme (RE) fragment level in order to leverage Hi-C data at its highest possible and native resolution. Briefly, HiCnv first computes 1D read coverage for each RE fragment, followed by normalization for GC content, mappability and fragment length, and by smoothing using kernel density estimation (KDE). KDE smoothed counts are divided into potential CNV segments using a Hidden Markov Model (HMM), and these segments are further processed for refinement of their breakpoint coordinates (segment ends) and assignment of their CNV labels. See our [paper](https://academic.oup.com/bioinformatics/article/34/2/338/4557186) for more details.

The current version of this HiCnv uses [Snakemake](https://snakemake.readthedocs.io/en/stable/) to facilitate deployment and improve reproducibility. Currently, HiCnv is setup to download SRA files followed by alignment with HiCPro and CNV analysis. PLEASE READ: Snakemake uses a system of templating which I will be using in this documentation, templating refers to the use of `{variable_name}` to generalize a file path and common templates I will use include:
- `{re}` - restriction enzyme
- `{srr}` - SRR ID
- `{cline}` - cell line 

## Authors

mod: jreyna@lji.org (Joaquin Reyna)

orig: abhijit@lji.org (Abhijit Chakraborty)

## Software dependencies
Rules within this workflow attempt to facilitate the installation of a few key software but there may still be some user-specific installation steps. The following is a list of required software and packages where applicable:
- Python 2.7 
- Python 3.7
  - snakemake 
  - ...
- R 
  - ...
- HiCPro - Git download the HiCPro repository using: `git clone https://github.com/nservant/HiC-Pro1
- HiCPro Docker img (installed by workflow)
- Singularity ([link](https://singularity.lbl.gov/install-linux))

## Usage
### Snakemake config
The Snakemake config file is located within `config/config.yaml`, please set the following variables (using the YAML format already provided) with paths to the correponding resources or software:
<pre>
R4: &#60;path to R 4.X.X with library dependencies installed&#62;
python2: &#60;path to Python 2.7 with package dependencies installed&#62;
hicpro_dir: &#60;path to the main HiCPro installation directory&#62;
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
and includes setting `GENOME_FRAGMENT = results/refs/restriction_enzymes/hg38_{re}_digestion.bed`. Once completed store the configuration file within `/results/refs/hicpro/config-hicpro.{re}.txt`. This process should be completed for every restriction enzyme used by your dataset.

### Download SRA paired fastq data
Uses Rule download_paired_fastq_sra:
<pre>
snakemake --profile workflow/profiles/pbs-torque results/main/{cline}/sra/{srr}_1.fastq.gz
</pre>

### Alignment
Process your Hi-C fastq files with [HiCPro pipeline](https://github.com/nservant/HiC-Pro) (Rule hicpro_align_only):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/main/{cline}/hicpro/bowtie_results/bwt2/{srr}/{srr}_1_hg38.bwt2merged.bam
</pre>

### Generate the coverage, GC content, mappability and fragment length information file
Download the mappability file which is used to generate the F_GC_MAP file (Rule download_hg38_mappability):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.sorted.bedGraph
</pre>

Rule process_refeature: generate a restriction fragment specific file known as the
*.fragments.F_GC_MAP.bed (Fragment length, GC content and Mappability information file
file using existings commands.
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
snakemake --profile workflow/profiles/pbs-torque results/main/{cline}/hicnv/{cline}_{srr}_hicnv_final.test
</pre>

### Run the HiCnv workflow the _fast way_
After installing all the necessary software, and setting up the configurations files you can simply run (Rule run_hicnv):
<pre>
snakemake --profile workflow/profiles/pbs-torque results/main/{cline}/hicnv/{cline}_{srr}_hicnv_final.test
</pre>
Snakemake deduces what rules need to be run and will run as if you you followed the slower tutorial from above.

## Running on a cluster system
If you are using PBS-Torque then you are in luck, this workflow has been set up to work to automatically submit each rule as a PBS-Torque job. To include your job submission system you must create a profile as specified [here](https://snakemake.readthedocs.io/en/stable/executing/cluster.html). To run with this workflow make sure to save your profile under `workflow/profiles`. Ideally, from the main directory you would do:
<pre>
cp -r /home/jreyna/jreyna/projects/HiCnv/workflow/profiles/pbs-torque /home/jreyna/jreyna/projects/HiCnv/workflow/profiles/<your new profile>
</pre>
Next, you would modify `/home/jreyna/jreyna/projects/HiCnv/workflow/profiles/<your new profile>/qsub.py` and more specificially you would change 
<pre>
cmd = 'qsub -l walltime=200:00:00,mem={}mb,nodes={}:ppn={} -o {} -e {} {}'
</pre>
so as to fit your job management system.

## Workflow schema
The following dag diagram was generated by snakemake and shows the connections between snakemake rules which for the basis for the HiCnv workflow.
<p align="center">
  <img src="https://raw.githubusercontent.com/joreynajr/HiCnv/master/images/graph.svg" width="800"/>
</p>
