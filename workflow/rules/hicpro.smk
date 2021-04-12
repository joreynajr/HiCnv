
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


# get the restriction enzyme digestion files
def re_digestion_file(wildcards):
    re = SAMPLESHEET.loc[wildcards.cline, 're']
    config = 'results/refs/restriction_enzymes/hg38_{}_digestion.bed'.format(re)
    return(config)


# get the restriction enzyme hicpro config
def re_config_file(wildcards):
    re = SAMPLESHEET.loc[wildcards.cline, 're']
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
        bowtie_results = directory('results/main/{cline}/hicpro/bowtie_results/bwt2/{cline}/'),
        bowtie_complete = touch('results/main/{cline}/hicpro/bowtie_results/{cline}/bowtie.complete')
    params:
        datadir1 = 'results/main/{cline}/hicpro/reads_syms/',
        datadir2 = 'results/main/{cline}/hicpro/reads_syms/{cline}/',
        outdir = 'results/main/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 80000
    log:
        'results/main/{cline}/logs/rule_hicpro_align_only_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_align_only_{cline}.bmk'
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


# Quality check the HiC data with HiCPro (single step)
rule hicpro_quality_checks_only: #  merging update incomplete
    input:
        bams = rules.hicpro_align_only.output.bowtie_complete,
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


## Process the HiC data with HiCPro (single step)
rule hicpro_hic_proc_only: # merging update complete
    input:
        bams = rules.hicpro_align_only.output.bowtie_complete,
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        vp = directory('results/main/{cline}/hicpro/hic_results/data/{cline}/'),
        vp_complete = touch('results/main/{cline}/hicpro/hic_results/data/{cline}/process.complete')
    params:
        datadir = 'results/main/{cline}/hicpro/bowtie_results/bwt2/',
        outdir = 'results/main/{cline}/hicpro/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 10000
    log:
        'results/main/{cline}/logs/rule_hicpro_hic_proc_only_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_hic_proc_only_{cline}.bmk'
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
        vp_complete = rules.hicpro_hic_proc_only.output.vp_complete,
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
# mat = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/10000/{cline}_10000.matrix',
# bed = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/10000/{cline}_10000_abs.bed'

rule hicpro_build_contact_maps_only: #  merging update complete
    input:
        merged_samples = rules.hicpro_merge_persample_only.output.vp,
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        mats = directory('results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/'),
        mats_complete = touch('results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/mats.complete')
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
#mat = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/10000/{cline}_10000.matrix',
#bed = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/raw/10000/{cline}_10000_abs.bed',
# matrix = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/iced/10000/{cline}_10000_iced.matrix',
# biases = 'results/main/{cline}/hicpro/hic_results/matrix/{cline}/iced/10000/{cline}_10000_iced.matrix.biases'
rule hicpro_ice_norm_only: #  merging update complete
    input:
        mats_complete = rules.hicpro_build_contact_maps_only.output.mats_complete,
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        iced = directory('results/main/{cline}/hicpro/hic_results/matrix/{cline}/iced/'),
        iced_complete = touch('results/main/{cline}/hicpro/hic_results/matrix/{cline}/iced/iced.complete')
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


## Run fastq_pair to ensure file is processed correctly
# Current template example can be found at:
# /mnt/BioHome/jreyna/archive/test_hicpro/test_data/dixon_2M/fastq_pair.sh
rule fastq_pair:
    input:
        unpack(get_r1_r2_fastqs),
    output:
        fastq_pair_dir = directory('results/main/{cline}/reads/fastq_pair/'),
        fastq_pair_flag = touch('results/main/{cline}/reads/fastq_pair/fastq_pair.complete')
    params:
        orig_fastq_dir = 'results/main/{cline}/reads/'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 64000
    log:
        'results/main/{cline}/logs/rule_fastq_pair_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_fastq_pair_{cline}.bmk'
    shell:
        """
            mkdir -p {output.fastq_pair_dir}

            # Uncompress r1 and r2 fastq's
            echo "# Uncompress r1 and r2 fastq's"
            for fn in {input.r1s} {input.r2s};
            do
                uncomp_fn=$(basename $fn | sed "s/\.gz//")
                uncomp_fn="{output.fastq_pair_dir}/$uncomp_fn"
                unpigz -p {resources.ppn} -k -c $fn > $uncomp_fn 2>> {log}
            done

            # run fastq_pair for all pairs
            echo "# run fastq_pair for all pairs"
            for fn in {input.r1s};
            do
                srr=$(echo $fn | xargs basename | cut -f 1 -d "_")
                fq1="{output.fastq_pair_dir}/${{srr}}_1.fastq"
                fq2="{output.fastq_pair_dir}/${{srr}}_2.fastq"
                fastq_pair $fq1 $fq2 >> {log} 2>&1
                rm $fq1 $fq2 >> {log} 2>&1
            done

            # compress the results
            echo "# compress the results"
            for fn in $(ls {output.fastq_pair_dir});
            do
                # compress
                rel_fn="{output.fastq_pair_dir}/${{fn}}"
                pigz -p {resources.ppn} $rel_fn 2>> {log}

                # rename the compressed file
                pigz_fn="{output.fastq_pair_dir}/${{fn}}.gz"
                final_fn=$(echo $fn | sed "s/fastq\.//" | sed "s/fq/fastq.gz/")
                final_fn="{output.fastq_pair_dir}/${{final_fn}}"
                mv $pigz_fn $final_fn 2>> {log}

                # remove uncompressed files
                #rm "{output.fastq_pair_dir}/${{fn}}" 2>> {log}
            done
        """


# Helper function to obtain all for a given cell line
# Make sure to download the accession list from SRA
def get_r1_r2_fastqs_with_adj(wildcards):

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
                r1 = 'results/main/{cline}/reads/fastq_pair/{acc}_1.paired.fastq.gz'.format(cline=wildcards.cline, acc=acc)
                r2 = 'results/main/{cline}/reads/fastq_pair/{acc}_2.paired.fastq.gz'.format(cline=wildcards.cline, acc=acc)
                r1s.append(r1)
                r2s.append(r2)
    d = {'r1s': r1s, 'r2s': r2s}
    return(d)


use rule hicpro_align_only as hicpro_align_only_with_adj with:
    input:
        unpack(get_r1_r2_fastqs_with_adj),
        fastq_pair_flag = 'results/main/{cline}/reads/fastq_pair/fastq_pair.complete',
        gs = rules.download_hg38_files.output.genome_sizes,
        digestion = re_digestion_file,
        bowtie2_idxs = rules.bowtie2_index_ref_genome.output,
        config = ancient(re_config_file),
        hicpro_img = rules.download_hicpro_singularity_img.output.hicpro_img
    output:
        bowtie_results = directory('results/main/{cline}/hicpro_with_adj/bowtie_results/bwt2/{cline}/'),
        bowtie_complete = touch('results/main/{cline}/hicpro_with_adj/bowtie_results/{cline}/bowtie.complete')
    params:
        datadir1 = 'results/main/{cline}/hicpro_with_adj/reads_syms/',
        datadir2 = 'results/main/{cline}/hicpro_with_adj/reads_syms/{cline}/',
        outdir = 'results/main/{cline}/hicpro_with_adj/'
    log:
        'results/main/{cline}/logs/rule_hicpro_align_only_with_adj_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_align_only_with_adj_{cline}.bmk'
