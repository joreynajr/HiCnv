# Align the HiC data
# https://github.com/nservant/HiC-Pro
# 1*.bwt2merged.bam and 2.bwt2merged.bam
# bam1 = 'data/{cline}/hicpro/{cline}.{srr}.1.bwt2merged.bam',
# bam2 = 'data/{cline}/hicpro/{cline}.{srr}.2.bwt2merged.bam'
rule hicpro_align_only:
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
                    HiC-Pro -s mapping \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Process the HiC data with HiCPro (single step)
rule hicpro_hic_proc_only:
    input:
        bam = rules.download_paired_fastq_sra.output.r1,
        config = 'refs/hicpro/config-hicpro.txt'
    output:
        proc_bam = ''
    params:
        datadir = '',
        outdir = ''
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 10000
    log:
        'logs/rule_hicpro_hic_proc_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            yes | singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s proc_hic \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Quality check the HiC data with HiCPro (single step)
rule hicpro_quality_checks_only:
    input:
        bam = rules.download_paired_fastq_sra.output.r1,
        config = 'refs/hicpro/config-hicpro.txt'
    output:
        proc_bam = ''
    params:
        datadir = '',
        outdir = ''
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 10000
    log:
        'logs/rule_hicpro_hic_proc_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            yes | singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s quality_checks \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Build contact maps for the HiC data with HiCPro (single step)
rule hicpro_build_contact_maps_only:
    input:
        bam = rules.download_paired_fastq_sra.output.r1,
        config = 'refs/hicpro/config-hicpro.txt'
    output:
        proc_bam = ''
    params:
        datadir = '',
        outdir = ''
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 10000
    log:
        'logs/rule_hicpro_hic_proc_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            yes | singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s build_contact_maps \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Build contact maps for the HiC data with HiCPro (single step)
rule hicpro_ice_norm_only:
    input:
        bam = rules.download_paired_fastq_sra.output.r1,
        config = 'refs/hicpro/config-hicpro.txt'
    output:
        proc_bam = ''
    params:
        datadir = '',
        outdir = ''
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 10000
    log:
        'logs/rule_hicpro_hic_proc_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            yes | singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s ice_norm \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """
