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
