# HiCnv version 3

HiCnv is a pipeline to call CNVs from Hi-C data and is now using Snakemake to facilitate the usage and interface of HiCnv.

# Steps: Process the HiC data; generate coverage, GC content, mappability and fragment length information file

1) Rule download_paired_fastq_sra: Download raw HiC data 
snakemake --profiles profile/pbs-torque `data/{cline}/sra/{srr}_1.fastq.gz`

2) Rule digest_reference_genome: digest hg38 using Mboi and HindIII
snakemake --profiles profile/pbs-torque `completed documentation coming`

3) Follow instructions at within the HiCPro documentation to generate an appropriate config file
Github link (https://github.com/nservant/HiC-Pro)

4) Rule hicpro_align: Process your Hi-C fastq files with HiCPro pipeline (https://github.com/nservant/HiC-Pro)
snakemake --profiles profile/pbs-torque `data/{cline}/{srr}/hicpro/{cline}.{srr}.ran.flag`

5) Rule download_hg38_mappability: download the mappability file which is used to generate the F_GC_MAP file
snakemake --profiles profile/pbs-torque `completed documentation coming`

6) Rule process_refeature: generate a restriction fragment specific file known as the *.fragments.F_GC_MAP.bed (Fragment length, GC content and Mappability information file) 
file using existings commands. 
snakemake --profiles profile/pbs-torque `completed documentation coming`

7) One dimensionalize your HiC data using pre-existing "scripts Read_coverage_generation/run_1DReadCoverage.pl" and create the *.perREfragStats file.
snakemake command tbd

8) Run hicnv_v2.R, check the usage for more details.

```bash
Usage: hicnv_v2.R [options]


Options:
        --refeature=REFEATURE
                A five column restriction enzyme cutsite bed file with GC content,
                mappability, fragment length infomation

                <chr>   <start> <end>   <GC frequency>  <mappability>   <fragment length>


        --coverage=COVERAGE
                A bedGraph file with read coverage signal (5th column).

                Alternatively, the .perREfragStats file


        --gccutoff=GCCUTOFF
                GC content cutoff. Anything below <gccutoff> will be removed [Default is 0.2].


        --mapcutoff=MAPCUTOFF
                Mappability cutoff. Anything below <mapcutoff> will be removed [Default is 0.5].


        --fragcutoff=FRAGCUTOFF
                Fragment length cutoff. Anything below <fragcutoff> will be removed [Default is 150].

                For Hi-C experiments with 4bp cut enzyme, this value is 150, for 6bp
                enzymes this value should be 1000.


        --refchrom=REFCHROM
                Name of the reference chromosome for CNV detection.

                If no name is provided then HiCnv proceed to estimate proportion of
                interaction count or PIC to decide a reference chromosome.


        --bandwidth=BANDWIDTH
                Genomic distance bandwidth with which Kernel smoothing will be performed [Default 1Mb].


        --hmmstate=HMMSTATE
                Number of HMM states to be searched in each chromosome [Default 10].


        --threshold=THRESHOLD
                Threshold value to define amplification and deletion with respect to
                normal region [Default is 0.2 i.e. deviation of 20% from mean normal
                value will be labeled as CNV].


        --prefix=PREFIX
                All the files and folders will be created with this name.


        --cpu=CPU
                Number of threads to be used [Default 1].


        -h, --help
                Show this help message and exit
```

# Other snakemake rules that are currently working and will be organized soon:
- Rule download_hg38_files
- Rule bowtie2_index_ref_genome

# Note:

CNV calling requires GC content, mappability and fragment length information of every RE
fragments. The *.F_GC_MAP.bed file contains all these information.  The file can be
created using F_GC_MAP.file.sh script available under "scripts/F_GC_MAP_Files/" folder.
To create the file, please run the script inside the "scripts/F_GC_MAP_Files/" folder.
Also, change the variables as per the Hi-C experiment.

For more details check "Rscript hicnv_v2.R --help".

# Covert an aligned Hi-C sam file into HiCnv usable format:

The script "samToHiCProFormat.pl" under "scripts/" folder takes an aligned HiC file in
sam format. This can be a merged.sam (both forward and reverse reads merged/paired
format) file or in single format file where forward and reverse reads are mapped into
separate files (e.g forward.sam and reverse.sam).

The script an be run like the following

perl samToHiCProFormat.pl -format paired -sam_file merged.sam -strand 0 -chr 3 -pos 4 -mapq 5 -read_len 6 -read_pos 10 -out_file test

or 

perl samToHiCProFormat.pl -format single -sam_file forward.sam,reverse.sam -strand 2 -chr 3 -pos 4 -mapq 5 -read_len 6 -read_pos 10 -out_file test

Check for the full details by "perl samToHiCProFormat.pl -help"

Script will create test.mate1.sam and test.mate2.sam file. Then add the header lines using samtools

samtools view -bT hg19.fa test.mate1.sam > mate1.bam

samtools view -bT hg19.fa test.mate2.sam > mate2.bam

Then in the run_1DReadCoverage.pl script of HiCnv, change the following variable to 

$hic_bwt2_folder_FWD = "mate1.bam";

$hic_bwt2_folder_REV = "mate2.bam";

This will enable HiCnv to read the bam files.

The "samToHiCProFormat_Example" folder contains example data to convert Hi-C sam files
into HiCnv format.

Note: In forward.sam and reverse,sam files the total number of reads should be equal
and assumed that each row of the two files should represent the same read.

# Double Minute (DM) and Homogeneously Staining Regions (HSR) scanning

Run ./scripts/dm_hsr.r to scan DM and HSR regions. 

For more details check "Rscript ./scripts/dm_hsr.r --help"

# Contact

mod: jreyna@lji.org (Joaquin Reyna)

orig: abhijit@lji.org (Abhijit Chakraborty)
