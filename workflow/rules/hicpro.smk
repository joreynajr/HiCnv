# HiCPro documentation https://github.com/nservant/HiC-Pro

# Align the HiC data
rule hicpro_align_only:
    input:
        r1 = rules.download_paired_fastq_sra.output.r1,
        r2 = rules.download_paired_fastq_sra.output.r2,
        config = 'refs/hicpro/config-hicpro.hindiii.txt'
    output:
        bam1 = 'data/{cline}/hicpro/bowtie_results/bwt2/{srr}/{srr}_1_hg38.bwt2merged.bam',
        bam2 = 'data/{cline}/hicpro/bowtie_results/bwt2/{srr}/{srr}_2_hg38.bwt2merged.bam'
    params:
        datadir1 = 'data/{cline}/hicpro/hicpro_tmp_{srr}/',
        datadir2 = 'data/{cline}/hicpro/hicpro_tmp_{srr}/{srr}/',
        outdir = 'data/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 64000
    log:
        'logs/rule_hicpro_align_only_{cline}_{srr}.log'
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
            singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s mapping \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1

            # removing the temp data dir
            # rm -r {params.datadir1}
        """


## Process the HiC data with HiCPro (single step)
rule hicpro_hic_proc_only:
    input:
        bam1 = rules.hicpro_align_only.output.bam1,
        bam2 = rules.hicpro_align_only.output.bam2,
        config = 'refs/hicpro/config-hicpro.hindiii.txt'
    output:
        vp = 'data/{cline}/hicpro/hic_results/data/{srr}/{srr}_hg38.bwt2pairs.validPairs'
    params:
        datadir = 'data/{cline}/hicpro/bowtie_results/bwt2/',
        outdir = 'data/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 10000
    log:
        'logs/rule_hicpro_hic_proc_only_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s proc_hic \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Quality check the HiC data with HiCPro (single step)
rule hicpro_quality_checks_only:
    input:
        bam1 = rules.hicpro_align_only.output.bam1,
        bam2 = rules.hicpro_align_only.output.bam2,
        config = 'refs/hicpro/config-hicpro.hindiii.txt'
    output:
        qc = directory('data/{cline}/hicpro/hic_results/pic/{srr}')
    params:
        datadir = 'data/{cline}/hicpro/bowtie_results/bwt2/',
        outdir = 'data/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 10000
    log:
        'logs/rule_hicpro_quality_checks_only_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s quality_checks \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Build contact maps for the HiC data with HiCPro (single step)
rule hicpro_merge_persample_only:
    input:
        vp = 'data/{cline}/hicpro/hic_results/data/{srr}/{srr}_hg38.bwt2pairs.validPairs',
        config = 'refs/hicpro/config-hicpro.hindiii.txt'
    output:
        vp = 'data/{cline}/hicpro/hic_results/data/{srr}/{srr}.allValidPairs',
        stats = 'data/{cline}/hicpro/hic_results/stats/{srr}/{srr}_allValidPairs.mergestat'
    params:
        datadir = 'data/{cline}/hicpro/hic_results/data/',
        outdir = 'data/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 64000
    log:
        'logs/rule_hicpro_merge_persample_only_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s merge_persample \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Build contact maps for the HiC data with HiCPro (single step)
rule hicpro_build_contact_maps_only:
    input:
        vp = 'data/{cline}/hicpro/hic_results/data/{srr}/{srr}.allValidPairs',
        config = 'refs/hicpro/config-hicpro.hindiii.txt'
    output:
        mat = 'data/{cline}/hicpro/hic_results/matrix/{srr}/raw/10000/{srr}_10000.matrix',
        bed = 'data/{cline}/hicpro/hic_results/matrix/{srr}/raw/10000/{srr}_10000_abs.bed'
    params:
        datadir = 'data/{cline}/hicpro/hic_results/data/',
        outdir = 'data/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 64000
    log:
        'logs/rule_hicpro_build_contact_maps_only_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s build_contact_maps \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Build contact maps for the HiC data with HiCPro (single step)
rule hicpro_ice_norm_only:
    input:
        mat = 'data/{cline}/hicpro/hic_results/matrix/{srr}/raw/10000/{srr}_10000.matrix',
        bed = 'data/{cline}/hicpro/hic_results/matrix/{srr}/raw/10000/{srr}_10000_abs.bed',
        config = 'refs/hicpro/config-hicpro.hindiii.txt'
    output:
        matrix = 'data/{cline}/hicpro/hic_results/matrix/{srr}/iced/10000/{srr}_10000_iced.matrix',
        biases = 'data/{cline}/hicpro/hic_results/matrix/{srr}/iced/10000/{srr}_10000_iced.matrix.biases'
    params:
        datadir = 'data/{cline}/hicpro/hic_results/matrix/',
        outdir = 'data/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 64000
    log:
        'logs/rule_hicpro_ice_norm_{cline}_{srr}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec software/hicpro_latest_ubuntu.img \
                    HiC-Pro -s ice_norm \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """
