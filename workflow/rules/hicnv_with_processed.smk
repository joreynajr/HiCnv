########## Running processed data like 4DN ##########

def get_4dn_bam_url(wildcards):
    # loading 4dn meta data
    meta = read_table('config/4dn_meta.tsv')
    meta = meta[meta.file_type == 'bam']
    indexes = []
    for i, sr in meta.iterrows():
        idx = '{}-B{}-T{}'.format(sr.exp_acc, sr.biorep, sr.techrep)
        indexes.append(idx)
    meta.index = indexes

    # getting the file accession for the current data
    file_acc = meta.loc[wildcards.srr, 'file_acc']
    url = 'https://data.4dnucleome.org/files-processed/{srr}/@@download/{srr}.bam'
    url = url.format(srr=file_acc)
    return(url)

# Download 4D Nucleome bam data
rule download_4d_nucleome_alns:
    params:
        key = 'GBYDTYDO',
        secret = 'kkrklx5k5n7v7bci',
        url = get_4dn_bam_url 
    output:
        'results/main/{cline}/4d_nucleome/{srr}.bam'
    log:
        'results/main/{cline}/logs/rule_download_4d_nucleome_alns_{cline}_{srr}.bam'
    shell:
        r"""
            curl -L \
                --user {params.key}:{params.secret} \
                -o {output} \
                {params.url} >> {log} 2>&1
        """


# Split an interleavened bam into R1 and R2 files
rule split_interleavened_bam:
    shell:
        r"""
            {config[samtools]} view -hbf 64 {input} > {output.r1} 2> {log}
            {config[samtools]} view -hbf 128 {input} > {output.r2} 2> {log}
        """

# Split the processed interleavened bam into R1 and R2 files
use rule split_interleavened_bam as split_processed_bam with:
    input:
        'results/main/{cline}/4d_nucleome/{srr}.bam'
    output:
        r1 = 'results/main/{cline}/4d_nucleome/{srr}_R1.bam',
        r2 = 'results/main/{cline}/4d_nucleome/{srr}_R2.bam'
    log:
        'results/main/{cline}/logs/rule_split_processed_bam_{cline}_{srr}.log'


# Split the processed interleavened bam into R1 and R2 files
rule sort_and_index_processed_bams:
    input:
        r1 = rules.split_processed_bam.output.r1,
        r2 = rules.split_processed_bam.output.r2
    output:
        r1 = 'results/main/{cline}/4d_nucleome/{srr}_R1.sorted.bam',
        r2 = 'results/main/{cline}/4d_nucleome/{srr}_R2.sorted.bam'
    log:
        'results/main/{cline}/logs/rule_split_processed_bam_{cline}_{srr}.log'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 50000
    shell:
        r"""
            # sort and index r1
            {config[samtools]} sort -o {output.r1} \
                                    -O BAM \
                                    -@ {resources.ppn} \
                                    {input.r1}
            {config[samtools]} index -@ {resources.ppn} {output.r1}

            # sort and index r2
            {config[samtools]} sort -o {output.r2} \
                                    -O BAM \
                                    -@ {resources.ppn} \
                                    {input.r2}
            {config[samtools]} index -@ {resources.ppn} {output.r2}
        """


# Make perREfragStats for processed data
rule make_perREfragStats_for_processed_alns_bioreps:
    input:
        fwd_alns = rules.sort_and_index_processed_bams.output.r1,
        rev_alns = rules.sort_and_index_processed_bams.output.r2,
        frag = re_digestion_file,
        fgc_map = re_fgc_map_file
    output:
        frag_stats = 'results/main/{cline}/hicnv/re_frag_stats/{cline}.{srr}.perREfragStats'
    params:
        merged_singletons = 'results/main/{cline}/hicnv/re_frag_stats/{cline}.{srr}.bwt2pairs.withSingles.mapq30.bam',
        outdir = 'results/main/{cline}/hicnv/re_frag_stats/{cline}_allMap2FragmentsOutput',
        prefix = '{cline}_allMap2FragmentsOutput'
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 50000
    log:
        'results/main/{cline}/logs/rule_make_perREfragStats_for_4dn_reps_{cline}_{srr}.log'
    shell:
        r"""
              # Merge singletons
              # pysam will through some warnings about a missing index
              # but this can be ignored according for forums
              {config[python2]} workflow/scripts/mergeSAM-singletons.py \
                            -f {input.fwd_alns} \
                            -r {input.rev_alns} \
                            -o {params.merged_singletons} \
                            -v \
                            --single \
                            -q 30 >> {log} 2>&1

              mkdir -p {params.outdir}

              # Map HiC fragments
              {config[python2]} workflow/scripts/mapped_2hic_fragments.py \
                            -f {input.frag} \
                            -r {params.merged_singletons} \
                            -s 100 \
                            -l 800 \
                            -d 1000 \
                            --all \
                            -v \
                            -o {params.outdir} >> {log} 2>&1

              # Sort the REfragStats
              cat {params.outdir}/{wildcards.cline}.{wildcards.srr}.bwt2pairs.withSingles.mapq30.perREfragStats \
                            | sort -k1,1 -k2,2n > {output} 2> {log}
          """


# Filtering for chromosomes 1-22 + X.
# The main hicnv_v2.R script will not work because there are too
# few datapoints generated for chrM and chrY should also be excluded.
# Saved in the combined tech_combined folder to continue running with the 
# rest of the pipeline. After this rules combine_bio_reps_perREfragStats,
# and run_bio_hicnv_auto_bandwidth will complete running HiCNV.
use rule filter_perREfragStats_for_main_chrs as filter_perREfragStats_for_processed_alns_bioreps with:
    input:
        frag_stats = rules.make_perREfragStats_for_processed_alns_bioreps.output.frag_stats,
        chr_list = 'resources/chromosome_lists/main_chrs.txt'
    output:
        frag_stats = 'results/main/{cline}/hicnv/tech_combined/{cline}.{srr}.perREfragStats'
    log:
        'results/main/{cline}/logs/rule_filter_perREfragStats_for_4dn_reps_{cline}_{srr}.log'
