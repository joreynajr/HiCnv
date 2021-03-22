configfile: 'config.yaml'
#report: "report/workflow.rst"

# Downloading hg38 files as needed
rule download_hg38_files:
    params:
        centromeres = 'refs/hg38/centromeres.txt.gz',
        gaps = 'refs/hg38/gap.txt.gz'
    output:
        genome = 'refs/hg38/hg38.fa.gz',
        genome_unzipped = 'refs/hg38/hg38.fa',
        genome_sizes = 'refs/hg38/hg38.chrom.sizes',
        centromeres = 'refs/hg38/centromeres.txt',
        gaps = 'refs/hg38/gap.txt'
    log:
        'logs/rule_download_hg38_files.out'
    shell:
        """
            wget -O {output.genome} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.fa.gz 2> {log}
            wget -O {output.genome_sizes} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/bigZips/hg38.chrom.sizes 2> {log}
            wget -O {params.centromeres} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/centromeres.txt.gz 2> {log}
            wget -O {params.gaps} https://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/gap.txt.gz 2> {log}

            # gunzip the centromeres and gaps
            gunzip --stdout {output.genome} > {output.genome_unzipped} 2> {log}
            gunzip {params.centromeres} 2> {log}
            gunzip {params.gaps} 2> {log}
        """


# Index the reference genome (GRCh38)
# http://bowtie-bio.sourceforge.net/bowtie2/manual.shtml#indexing-a-reference-genome
rule bowtie2_index_ref_genome:
    input:
        rules.download_hg38_files.output.genome_unzipped
    output:
        multiext('refs/hg38/hg38', '.1.bt2', '.2.bt2', '.3.bt2', '.4.bt2', '.rev.1.bt2', '.rev.2.bt2')
    params:
        base_out = 'refs/hg38/hg38'
    log:
        'logs/rule_bowtie2_index_ref_genome.log'
    shell:
        """
            bowtie2-build {input} {params.base_out} >> {log} 2>&1
        """


# Digesting the reference genome with different RE's, for now
# I am including MboI, and HindIII. I followed the following instructions:
# https://github.com/nservant/HiC-Pro/blob/master/doc/UTILS.md.
rule digest_reference_genome:
    input:
        'refs/hg38/hg38.fa'
    output:
        mboi = 'refs/restriction_enzymes/hg38_mboi_digestion.bed',
        hindiii = 'refs/restriction_enzymes/hg38_hindiii_digestion.bed'
    log:
        'logs/rule_digest_reference_genome.log'
    shell:
        """
            {config[python2]} {config[hicpro_utils]}/digest_genome.py -r mboi -o {output.mboi} {input}
            {config[python2]} {config[hicpro_utils]}/digest_genome.py -r hindiii -o {output.hindiii} {input}
        """


# rule to download sra files which derived from paired sequencing
# using the grabseqs tool and modified code from the hisss snakemake pipeline
# https://github.com/louiejtaylor/hisss/blob/master/rules/sra_paired.rules
rule download_paired_fastq_sra:
    output:
        r1 = 'data/{cline}/sra/{srr}_1.fastq.gz',
        r2 = 'data/{cline}/sra/{srr}_2.fastq.gz',
        meta = 'data/{cline}/sra/{srr}.meta.csv'
    log:
        'logs/rule_download_paired_fastq_sra_{cline}_{srr}.log'
    shadow: 'shallow'
    params:
        outdir = 'data/{cline}/sra/',
        meta = '{srr}.meta.csv'
    resources:
        mem_mb = 8000,
        nodes = 1,
        ppn = 4
    shell:
        """
            grabseqs sra -m {params.meta} \
                        -o {params.outdir} \
                        -r 4 \
                        -t {resources.ppn} \
                        {wildcards.srr} >> {log} 2>&1
        """


# Align the HiC data
# https://github.com/nservant/HiC-Pro
# 1*.bwt2merged.bam and 2.bwt2merged.bam
# bam1 = 'data/{cline}/hicpro/{cline}.{srr}.1.bwt2merged.bam',
# bam2 = 'data/{cline}/hicpro/{cline}.{srr}.2.bwt2merged.bam'
rule hicpro_align:
    input:
        r1 = rules.download_paired_fastq_sra.output.r1,
        r2 = rules.download_paired_fastq_sra.output.r2,
        config = 'refs/hicpro/config-hicpro.txt'
    output:
        flag = 'data/{cline}/{srr}/hicpro/{cline}.{srr}.ran.flag'
    params:
        datadir1 = 'data/{cline}/hicpro_tmp_{srr}/',
        datadir2 = 'data/{cline}/hicpro_tmp_{srr}/{srr}/',
        outdir = 'data/{cline}/{srr}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 50000
    log:
        'logs/rule_hicpro_align_{cline}_{srr}.log'
    shell:
        """
            # setting up a temporary data directory structure for hicpro
            mkdir -p {params.datadir2}
            abs_r1=$(readlink -f {input.r1})
            abs_r2=$(readlink -f {input.r2})
            ln -s $abs_r1 {params.datadir2}
            ln -s $abs_r2 {params.datadir2}

            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir1})
            abs_outdir=$(readlink -f {params.outdir})

            # running without setting -s so that it runs the entire
            # pipeline (default settings)
            yes | singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1

            # removing the temp data dir
            #rm -r {params.datadir1}

            touch {output.flag}
        """

include: "rules/hicpro.smk"
include: "rules/hicnv.smk"
