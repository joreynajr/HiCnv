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
        'results/main/{cline}/4d_nucleome/{srr}/{srr}.bam'
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
# For more information check out : https://broadinstitute.github.io/picard/explain-flags.html
rule split_interleavened_bam:
    resources:
        nodes = 1,
        ppn = 2
    shell:
        r"""
            # parallelizing by giving each call a process
            {config[samtools]} view -hbf {params.flag1} {input} > {output.r1} 2> {log} & 
            {config[samtools]} view -hbf {params.flag2} {input} > {output.r2} 2> {log} & 
            wait 
        """


# Split the processed interleavened bam into R1 and R2 files
# flag1 = 64 means "first in pair" (only)
# flag2 = 128 means "second in pair" (only)
# For more information check out : https://broadinstitute.github.io/picard/explain-flags.html
use rule split_interleavened_bam as split_processed_bam_64_128 with:
    input:
        'results/main/{cline}/4d_nucleome/{srr}.bam'
    output:
        r1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R1.bam',
        r2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R2.bam'
    params:
        flag1 = 64,
        flag2 = 128
    log:
        'results/main/{cline}/logs/rule_split_processed_bam_{cline}_{srr}.log'


## Split the processed interleavened bam into R1 and R2 files
## flag1 = 65 means "read paired" & "first in pair" (only)
## flag2 = 129 means "read paired" & "second in pair" (only)
## This is a test to find out if the 65 and 129 flags can ensure 
## every n reads in R1 has a corresponding read in R2 (also n reads).
## For more information check out : https://broadinstitute.github.io/picard/explain-flags.html
#use rule split_interleavened_bam as split_processed_bam_65_129 with:
#    input:
#        'results/main/{cline}/4d_nucleome/{srr}.bam'
#    output:
#        r1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_65_129_R1.bam',
#        r2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_65_129_R2.bam'
#    params:
#        flag1 = 65,
#        flag2 = 129
#    log:
#        'results/main/{cline}/logs/rule_split_processed_bam_65_129_{cline}_{srr}.log'


## Sort and index the processed bams 
#rule sort_and_index_processed_bams:
#    input:
#        r1 = rules.split_processed_bam_64_128.output.r1,
#        r2 = rules.split_processed_bam_64_128.output.r2
#    output:
#        r1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R1.sorted.bam',
#        r2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R2.sorted.bam',
#        bai1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R1.sorted.bam.bai',
#        bai2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R2.sorted.bam.bai'
#    log:
#        'results/main/{cline}/logs/rule_split_processed_bam_{cline}_{srr}.log'
#    resources:
#        nodes = 1,
#        ppn = 4,
#        mem_mb = 50000
#    shell:
#        r"""
#            # sort and index r1
#            {config[samtools]} sort -o {output.r1} \
#                                    -O BAM \
#                                    -@ {resources.ppn} \
#                                    {input.r1}
#            {config[samtools]} index -@ {resources.ppn} {output.r1}
#
#            # sort and index r2
#            {config[samtools]} sort -o {output.r2} \
#                                    -O BAM \
#                                    -@ {resources.ppn} \
#                                    {input.r2}
#            {config[samtools]} index -@ {resources.ppn} {output.r2}
#        """

# Sort and index the processed bams 
rule query_sort_processed_bams:
    input:
        r1 = rules.split_processed_bam_64_128.output.r1,
        r2 = rules.split_processed_bam_64_128.output.r2
    output:
        r1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R1.query_sorted.bam',
        r2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R2.query_sorted.bam',
    log:
        'results/main/{cline}/logs/rule_query_sort_proccessed_bam_{cline}_{srr}.log'
    resources:
        nodes = 1,
        ppn = 8,
        mem_mb = 50000
    shell:
        r"""
            # query sort R1
            {config[samtools]} sort -o {output.r1} \
                                    -O BAM \
                                    -n \
                                    -@ {resources.ppn} \
                                    {input.r1} 2> {log}

            # query sort R2
            {config[samtools]} sort -o {output.r2} \
                                    -O BAM \
                                    -n \
                                    -@ {resources.ppn} \
                                    {input.r2} 2> {log}
        """

# Lineup the sorted R1 and R2 files
# Lineup means we find the query names that are the same 
# between R1 and R2 and extract only
rule lineup_processed_bams:
    input:
        r1 = rules.query_sort_processed_bams.output.r1,
        r2 = rules.query_sort_processed_bams.output.r2
    output:
        r1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R1.sorted.linedup.bam',
        r2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R2.sorted.linedup.bam'
    log:
        'results/main/{cline}/logs/rule_lineup_processed_bams_{cline}_{srr}.log'
    resources:
        nodes = 1,
        ppn = 4,
        mem_mb = 4000
    shell:
        r"""
            # Lineup the sorted R1 and R2 files
            #{config[python_hicpro]} workflow/scripts/lineup_paired_reads.py \
            #{config[python_hicpro]} workflow/scripts/lineup_paired_reads_with_tries.py \
            {config[python_hicpro]} workflow/scripts/lineup_paired_reads_low_mem.py \
                                        --r1 {input.r1} \
                                        --r2 {input.r2} \
                                        --r1-out {output.r1} \
                                        --r2-out {output.r2} >> {log} 2>&1
        """


# non-sorted version
#r1 = rules.split_processed_bam.output.r1,
#r2 = rules.split_processed_bam.output.r2,

# sorted version 
#r1 = rules.sort_and_index_processed_bams.output.r1,
#r2 = rules.sort_and_index_processed_bams.output.r2,

# linedup version 
#r1 = rules.lineup_processed_bams.output.r1,
#r2 = rules.lineup_processed_bams.output.r2,

# Make perREfragStats for processed data
rule make_perREfragStats_for_processed_alns_bioreps:
    input:
        r1 = rules.lineup_processed_bams.output.r1,
        r2 = rules.lineup_processed_bams.output.r2,
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
                            -f {input.r1} \
                            -r {input.r2} \
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
