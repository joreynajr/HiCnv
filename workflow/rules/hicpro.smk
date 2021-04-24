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


# Renaming because _1 and _2 in file names can caused a problem
# which meant I reverted to using _R1 and _R2 in the HiC-Pro configuration file.
rule rename_before_hicpro:
    input:
        unpack(get_r1_r2_fastqs)
    output:
        new_dir = directory('results/main/{cline}/hicpro/renamed_fastqs/{cline}/'),
        rename_complete = touch('results/main/{cline}/hicpro/renamed_fastqs/renamed.complete')
    log:
        'results/main/{cline}/logs/rule_rename_before_hicpro_{cline}.log'
    shell:
        """
            mkdir -p {output.new_dir}

            # renaming R1's
            for fn in {input.r1s};
            do
                new_fn=$(basename $fn | sed "s/_1\.fastq\.gz/_R1.fastq.gz/")
                new_fn="{output.new_dir}/$new_fn"
                full_orig=$(readlink -f $fn)
                ln -s $full_orig $new_fn
            done

            # renaming R2's
            for fn in {input.r2s};
            do
                new_fn=$(basename $fn | sed "s/_2\.fastq\.gz/_R2.fastq.gz/")
                new_fn="{output.new_dir}/$new_fn"
                full_orig=$(readlink -f $fn)
                ln -s $full_orig $new_fn
            done
        """


# Align the HiC data with merging capability
rule hicpro_align_only: # merging update complete
    input:
        fastq_dir = rules.rename_before_hicpro.output.new_dir,
        gs = rules.download_hg38_files.output.genome_sizes,
        digestion = re_digestion_file,
        bowtie2_idxs = rules.bowtie2_index_ref_genome.output,
        config = ancient(re_config_file),
    output:
        bowtie_complete = touch('results/main/{cline}/hicpro/bowtie_results/bowtie.complete')
    params:
        datadir = 'results/main/{cline}/hicpro/renamed_fastqs/', # part of rule rename_before_hicpr
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
            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {params.datadir})
            abs_outdir=$(readlink -f {params.outdir})

            # running with setting -s so that it runs only the alignment 
            # pipeline (default settings)
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
    output:
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
                    HiC-Pro -s ice_norm \
                            -i $abs_datadir \
                            -o $abs_outdir \
                            -c {input.config} >> {log} 2>&1
        """
