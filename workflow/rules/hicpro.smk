import glob
from pandas import read_table

# HiCPro documentation https://github.com/nservant/HiC-Pro
# Download the HiCPro singularity img file which as linked
# under the HiCPro documentation section called
# "Using HiC-Pro through Singularity", it's the first sentence.
rule download_hicpro_singularity_img:
    output:
        hicpro_img = 'resources/software/hicpro_latest_ubuntu.img'
    shell:
        """
            wget -O {output} https://zerkalo.curie.fr/partage/HiC-Pro/singularity_images/hicpro_latest_ubuntu.img
        """


# load sample to re dict
def load_samplesheet():
    ss = read_table('config/samplesheet.tsv', index_col=0)
    return(ss)


# get the restriction enzyme digestion files
def re_digestion_file(wildcards):
    ss = load_samplesheet()
    re = ss.loc[wildcards.cline, 're']
    config = 'results/refs/restriction_enzymes/hg38_{}_digestion.bed'.format(re)
    return(config)


# get the restriction enzyme hicpro config
def re_config_file(wildcards):
    ss = load_samplesheet()
    re = ss.loc[wildcards.cline, 're']
    config = 'results/refs/hicpro/config-hicpro.{}.txt'.format(re)
    return(config)


# Helper function to obtain all for a given cell line
# Make sure to download the accession list from SRA
def get_r1_r2_fastqs(wildcards):

    # init lists to collect r1 and r2 files
    r1s = []
    r2s = []

    # list the accession files for this sample
    acc_lists = glob.glob('results/main/{cline}/reads/{cline}.*.SRR_Acc_List.txt'.format(cline=wildcards.cline))

    # parse through the accession lists and get teh r1 and r2 paths 
    for acc_list in acc_lists:
        with open(acc_list) as fr:
            accs = [x.strip() for x in fr.readlines()]
            for acc in accs:
                r1 = 'results/main/{cline}/reads/{acc}_1.fastq.gz'.format(cline=wildcards.cline, acc=acc)
                r2 = 'results/main/{cline}/reads/{acc}_2.fastq.gz'.format(cline=wildcards.cline, acc=acc)
                r1s.append(r1)
                r2s.append(r2)
    d = {'r1s': r1s, 'r2s': r2s}
    return(d)


# Helper function to obtain all for a given cell line
# Make sure to download the accession list from SRA
def get_bam1s_bam2s(wildcards):

    # init lists to collect bam1 and bam2 files
    bam1s = []
    bam2s = []

    # list the accession files for this sample
    acc_lists = glob.glob('results/main/{cline}/reads/{cline}.*.SRR_Acc_List.txt'.format(cline=wildcards.cline))

    # parse through the accession lists and get teh bam1 and bam2 paths
    for acc_list in acc_lists:
        with open(acc_list) as fr:
            accs = [x.strip() for x in fr.readlines()]
            for acc in accs:
                bam1 = 'results/main/{cline}/hicpro/bowtie_results/bwt2/{acc}_1_hg38.bwt2merged.bam'.format(cline=wildcards.cline, acc=acc)
                bam2 = 'results/main/{cline}/hicpro/bowtie_results/bwt2/{acc}_2_hg38.bwt2merged.bam'.format(cline=wildcards.cline, acc=acc)
                bam1s.append(bam1)
                bam2s.append(bam2)
    d = {'bam1s': bam1s, 'bam2s': bam2s}
    return(d)


# Helper function to obtain all for a given cell line
# # Make sure to download the sra accession list from SRA
def get_srrs(wildcards):
    # list the sra accession files for this sample
    srr_files = glob.glob('results/main/{cline}/reads/{cline}.*.SRR_Acc_List.txt'.format(cline=wildcards.cline))
    srrs = []

    # parse through the sra accession lists and get teh r1 and r2 paths
    for srr_file in srr_files:
        with open(srr_file) as fr:
            curr = [x.strip() for x in fr.readlines()]
            srrs.extend(curr)
    print(srrs)
    return(srrs)


# Align the HiC data with merging capability
rule hicpro_align_only: # merging update complete
    input:
        unpack(get_r1_r2_fastqs),
        gs = rules.download_hg38_files.output.genome_sizes,
        digestion = re_digestion_file,
        bowtie2_idxs = rules.bowtie2_index_ref_genome.output,
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        bam1 = 'results/main/{cline}/hicpro/bowtie_results/bwt2/{cline}/{srr}_1_hg38.bwt2merged.bam',
        bam2 = 'results/main/{cline}/hicpro/bowtie_results/bwt2/{cline}/{srr}_2_hg38.bwt2merged.bam'
    params:
        datadir1 = 'results/main/{cline}/hicpro/reads_syms/',
        datadir2 = 'results/main/{cline}/hicpro/reads_syms/{cline}/',
        outdir = 'results/main/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 80000
    log:
        'results/main/{cline}/logs/rule_hicpro_align_only_{cline}_{srr}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_align_only_{cline}_{srr}.bmk'
    shell:
        """
            # setting up a temporary data directory structure for hicpro
            mkdir -p {params.datadir2}

            for fn in {input.r1s};
            do
                abs_r1=$(readlink -f $fn)
                ln -f -s $abs_r1 {params.datadir2}
            done

            for fn in {input.r2s};
            do
                abs_r2=$(readlink -f $fn)
                ln -f -s $abs_r2 {params.datadir2}
            done

            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir1})
            abs_outdir=$(readlink -f {params.outdir})

            # running without setting -s so that it runs the entire
            # pipeline (default settings)
            singularity exec {input.hicpro_img} \
                    HiC-Pro -s mapping \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


## Align the HiC data (deprecated, was no merging samples)
#rule hicpro_align_only:
#    input:
#        r1 = rules.download_paired_fastq_sra.output.r1,
#        r2 = rules.download_paired_fastq_sra.output.r2,
#        gs = rules.download_hg38_files.output.genome_sizes,
#        digestion = re_digestion_file,
#        bowtie2_idxs = rules.bowtie2_index_ref_genome.output,
#        config = ancient(re_config_file),
#        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
#    output:
#        bam1 = 'results/main/{cline}/hicpro/bowtie_results/bwt2/{srr}/{srr}_1_hg38.bwt2merged.bam',
#        bam2 = 'results/main/{cline}/hicpro/bowtie_results/bwt2/{srr}/{srr}_2_hg38.bwt2merged.bam'
#    params:
#        datadir1 = 'results/main/{cline}/hicpro/hicpro_tmp_{srr}/',
#        datadir2 = 'results/main/{cline}/hicpro/hicpro_tmp_{srr}/{srr}/',
#        outdir = 'results/main/{cline}/hicpro/'
#    resources:
#        nodes = 1,
#        ppn = 4,
#        mem_mb = 80000
#    log:
#        'results/main/{cline}/logs/rule_hicpro_align_only_{cline}_{srr}.log'
#    shell:
#        """
#            # setting up a temporary data directory structure for hicpro
#            mkdir -p {params.datadir2}
#            abs_r1=$(readlink -f {input.r1})
#            abs_r2=$(readlink -f {input.r2})
#            ln -f -s $abs_r1 {params.datadir2}
#            ln -f -s $abs_r2 {params.datadir2}
#
#            # getting absoluate paths for data and outdirs, required
#            # by HiCPro
#            abs_datadir=$(readlink -f {params.datadir1})
#            abs_outdir=$(readlink -f {params.outdir})
#
#            # running without setting -s so that it runs the entire
#            # pipeline (default settings)
#            singularity exec {input.hicpro_img} \
#                    HiC-Pro -s mapping \
#                            -i $abs_datadir \
#                            -o $abs_outdir \
#                            -c {input.config} >> {log} 2>&1
#
#            # removing the temp data dir
#            # rm -r {params.datadir1}
#        """



## Process the HiC data with HiCPro (single step)
rule hicpro_hic_proc_only: # merging update complete
    input:
        bams = rules.hicpro_align_only.output,
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        vp = 'results/main/{cline}/hicpro/hic_results/data/{cline}/{srr}_hg38.bwt2pairs.validPairs'
    params:
        datadir = 'results/main/{cline}/hicpro/bowtie_results/bwt2/',
        outdir = 'results/main/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 10000
    log:
        'results/main/{cline}/logs/rule_hicpro_hic_proc_only_{cline}_{srr}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_hic_proc_only_{cline}_{srr}.bmk'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec {input.hicpro_img} \
                    HiC-Pro -s proc_hic \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Quality check the HiC data with HiCPro (single step)
rule hicpro_quality_checks_only: #  merging update incomplete
    input:
        bams = rules.hicpro_align_only.output,
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        qc = directory('results/main/{cline}/hicpro/hic_results/pic/{cline}')
    params:
        datadir = 'results/main/{cline}/hicpro/bowtie_results/bwt2/',
        outdir = 'results/main/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 10000
    log:
        'results/main/{cline}/logs/rule_hicpro_quality_checks_only_{cline}.log'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec {input.hicpro_img} \
                    HiC-Pro -s quality_checks \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Helper function to obtain all for a given cell line
# Make sure to download the accession list from SRA
def get_vps(wildcards):

    # init lists to collect bam1 and bam2 files
    vps = []

    # list the accession files for this sample
    acc_lists = glob.glob('results/main/{cline}/reads/{cline}.*.SRR_Acc_List.txt'.format(cline=wildcards.cline))

    # parse through the accession lists and get teh bam1 and bam2 paths
    for acc_list in acc_lists:
        with open(acc_list) as fr:
            accs = [x.strip() for x in fr.readlines()]
            for acc in accs:
                vp = 'results/main/{cline}/hicpro/hic_results/data/{cline}/{acc}_hg38.bwt2pairs.validPairs'.format(cline=wildcards.cline, acc=acc)
                vps.append(vp)
    return(vps)


# Build contact maps for the HiC data with HiCPro (single step)
rule hicpro_merge_persample_only: #  merging update complete
    input:
        vp = get_vps,
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        vp = 'results/main/{cline}/hicpro/hic_results/data/{cline}/{cline}.allValidPairs',
        stats = 'results/main/{cline}/hicpro/hic_results/stats/{cline}/{cline}_allValidPairs.mergestat'
    params:
        datadir = 'results/main/{cline}/hicpro/hic_results/data/',
        outdir = 'results/main/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 64000
    log:
        'results/main/{cline}/logs/rule_hicpro_merge_persample_only_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_merge_persample_only_{cline}.bmk'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec {input.hicpro_img} \
                    HiC-Pro -s merge_persample \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Build contact maps for the HiC data with HiCPro (single step)
rule hicpro_build_contact_maps_only: #  merging update complete
    input:
        vp = 'results/main/{cline}/hicpro/hic_results/data/{cline}/{cline}.allValidPairs',
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        mat = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/10000/{cline}_10000.matrix',
        bed = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/10000/{cline}_10000_abs.bed'
    params:
        datadir = 'results/main/{cline}/hicpro/hic_results/data/',
        outdir = 'results/main/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 64000
    log:
        'results/main/{cline}/logs/rule_hicpro_build_contact_maps_only_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_build_contact_maps_only_{cline}.bmk'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec {input.hicpro_img} \
                    HiC-Pro -s build_contact_maps \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """


# Build contact maps for the HiC data with HiCPro (single step)
rule hicpro_ice_norm_only: #  merging update complete
    input:
        mat = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/10000/{cline}_10000.matrix',
        bed = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/10000/{cline}_10000_abs.bed',
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        matrix = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/iced/10000/{cline}_10000_iced.matrix',
        biases = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/iced/10000/{cline}_10000_iced.matrix.biases'
    params:
        datadir = 'results/main/{cline}/hicpro/hic_results/matrix/',
        outdir = 'results/main/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 64000
    log:
        'results/main/{cline}/logs/rule_hicpro_ice_norm_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_ice_norm_{cline}.bmk'
    shell:
        """
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # pipeline (default settings)
            singularity exec {input.hicpro_img} \
                    HiC-Pro -s ice_norm \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """
