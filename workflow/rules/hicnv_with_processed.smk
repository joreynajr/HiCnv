#########################################################################################
#########################################################################################
############################ Running processed data like 4DN ############################
#########################################################################################
#########################################################################################

# Helper function for one dimensionalize_pairs
def get_digestion_file(wildcards):
    re = SAMPLESHEET.loc[wildcards.cline, 're']
    config = 'results/refs/restriction_enzymes/hg38_{}_digestion.bed'.format(re)
    return(config)

rule one_dimensionalize_pairs_into_coverage_bed: # (status: running)
    input:
        pairs = rules.download_4d_nucleome_pairs.output,
        re = get_digestion_file
    params:
        prefix =  'results/main/{cline}/4d_nucleome/{srr}/{srr}',
        bedtools = '/mnt/BioApps/bedtools/bin/'
    output:
        coverage = protected('results/main/{cline}/4d_nucleome/tech_reps/{srr}/{srr}.coverage.sorted.bed')
    log:
        'results/main/{cline}/logs/rule_one_dimensionalize_pairs_into_coverage_bed_{cline}_{srr}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_one_dimensionalize_pairs_into_coverage_bed_{cline}_{srr}.bmk'
    resources:
        mem_mb = 100000
    shell:
        r"""
            {config[python_hicpro]} workflow/scripts/one_dimensionalize_pairs.py \
                    --pairs {input.pairs} \
                    --re {input.re} \
                    --prefix {params.prefix} \
                    --bedtools {params.bedtools} >> {log} 2>&1
        """


rule exclude_non_hicnv_chroms_from_coverage: # (status: running)
    input:
        rules.one_dimensionalize_pairs_into_coverage_bed.output.coverage
    output:
        'results/main/{cline}/4d_nucleome/tech_reps/{srr}/{srr}.coverage.sorted.fltchr.bed'
    log:
        'results/main/{cline}/logs/rule_exclude_non_hicnv_chroms_from_coverage_{cline}_{srr}.log'
    shell:
        r"""
            # use perl regex to extract normal chromosomes (1-22 + X)
            # potential bug if user introduces new chromsomes with new numbers
            grep -P 'chr[0-9X]{{1,2}}(?=\s)' {input} > {output}
        """


rule check_final_coverage_and_re: #(status: running)
    input:
        coverage = rules.exclude_non_hicnv_chroms_from_coverage.output,
        re = get_digestion_file
    output:
        check = 'results/main/{cline}/4d_nucleome/tech_reps/{srr}/{srr}.coverage.check',
        diff = 'results/main/{cline}/4d_nucleome/tech_reps/{srr}/{srr}.coverage.diff'
    log:
        'results/main/{cline}/logs/rule_check_final_coverage_and_re_{cline}_{srr}.log'
    shell:
        r"""
            # check 1: count the number of lines in both files
            echo "# check 1: count the number of lines in both files" >> {log}
            echo 'wc -l for the coverage file' >> {output.check}
            cov_check="{input.coverage}.input.check"
            cat {input.coverage} | cut -f 1-3 > $cov_check
            cat $cov_check | wc -l >> {output.check}

            echo >> {output.check}
            echo 'wc -l for the re file' >> {output.check}
            re_check="{input.re}.input.check"
            cat {input.re} | grep -P 'chr[0-9X]{{1,2}}(?=\s)' | cut -f 1-3 > $re_check
            cat $re_check | wc -l >> {output.check}

            # check 2: diff both chr filtered files
            echo "# check 2: diff both chr filtered files" >> {log}
            diff $cov_check $re_check > {output.diff} && \
                echo "no difference" > {log} || echo "differences!" >> {log}
        """


# Warning: You will have a bug with this file structure, I recommened making a 
# technical replicates folder
def get_replicate_coverage_fns(wildcards):
    fns = 'results/main/{cline}/4d_nucleome/tech_reps/*/*.coverage.sorted.fltchr.bed'
    fns = fns.format(**wildcards)
    fns = glob.glob(fns)

    # return None if empty
    if len(fns) == 0:
        #raise MissingInputException(rules.map_replicates_4dn, 'results/main/{cline}/4d_nucleome/tech_reps/*/*.coverage.sorted.fltchr.bed')
        raise Exception('Missing results/main/{cline}/4d_nucleome/tech_reps/*/*.coverage.sorted.fltchr.bed')
    return(fns)

rule map_replicates_4dn: # (status: running)
    input:
        get_replicate_coverage_fns
    output:
        protected('results/main/{cline}/4d_nucleome/merged/{cline}.coverage.sorted.fltchr.bedgraph')
    log:
        'results/main/{cline}/logs/rule_map_replicates_4dn_{cline}.log'
    shell:
        r"""
            {config[python_hicpro]} workflow/scripts/map_multiple_identical_bedgraphs.py \
                       --bedgraphs {input} \
                       --out {output} \
                       --check
        """


## Split an interleavened bam into R1 and R2 files
## For more information check out : https://broadinstitute.github.io/picard/explain-flags.html
#rule split_interleavened_bam:
#    resources:
#        nodes = 1,
#        ppn = 2
#    shell:
#        r"""
#            # parallelizing by giving each call a process
#            {config[samtools]} view -hbf {params.flag1} {input} > {output.r1} 2> {log} & 
#            {config[samtools]} view -hbf {params.flag2} {input} > {output.r2} 2> {log} & 
#            wait 
#        """
#
#
## Split the processed interleavened bam into R1 and R2 files
## flag1 = 64 means "first in pair" (only)
## flag2 = 128 means "second in pair" (only)
## For more information check out : https://broadinstitute.github.io/picard/explain-flags.html
#use rule split_interleavened_bam as split_processed_bam_64_128 with:
#    input:
#        'results/main/{cline}/4d_nucleome/{srr}.bam'
#    output:
#        r1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R1.bam',
#        r2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R2.bam'
#    params:
#        flag1 = 64,
#        flag2 = 128
#    log:
#        'results/main/{cline}/logs/rule_split_processed_bam_{cline}_{srr}.log'
#
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
#
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
#
## Sort and index the processed bams 
#rule query_sort_processed_bams:
#    input:
#        r1 = rules.split_processed_bam_64_128.output.r1,
#        r2 = rules.split_processed_bam_64_128.output.r2
#    output:
#        r1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R1.query_sorted.bam',
#        r2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R2.query_sorted.bam',
#    log:
#        'results/main/{cline}/logs/rule_query_sort_proccessed_bam_{cline}_{srr}.log'
#    resources:
#        nodes = 1,
#        ppn = 8,
#        mem_mb = 50000
#    shell:
#        r"""
#            # query sort R1
#            {config[samtools]} sort -o {output.r1} \
#                                    -O BAM \
#                                    -n \
#                                    -@ {resources.ppn} \
#                                    {input.r1} 2> {log}
#
#            # query sort R2
#            {config[samtools]} sort -o {output.r2} \
#                                    -O BAM \
#                                    -n \
#                                    -@ {resources.ppn} \
#                                    {input.r2} 2> {log}
#        """
#
## Lineup the sorted R1 and R2 files
## Lineup means we find the query names that are the same 
## between R1 and R2 and extract only
#rule lineup_processed_bams:
#    input:
#        r1 = rules.query_sort_processed_bams.output.r1,
#        r2 = rules.query_sort_processed_bams.output.r2
#    output:
#        r1 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R1.sorted.linedup.bam',
#        r2 = 'results/main/{cline}/4d_nucleome/{srr}/{srr}_R2.sorted.linedup.bam'
#    log:
#        'results/main/{cline}/logs/rule_lineup_processed_bams_{cline}_{srr}.log'
#    resources:
#        nodes = 1,
#        ppn = 4,
#        mem_mb = 4000
#    shell:
#        r"""
#            # Lineup the sorted R1 and R2 files
#            #{config[python_hicpro]} workflow/scripts/lineup_paired_reads.py \
#            #{config[python_hicpro]} workflow/scripts/lineup_paired_reads_with_tries.py \
#            {config[python_hicpro]} workflow/scripts/lineup_paired_reads_low_mem.py \
#                                        --r1 {input.r1} \
#                                        --r2 {input.r2} \
#                                        --r1-out {output.r1} \
#                                        --r2-out {output.r2} >> {log} 2>&1
#        """
#
# non-sorted version
#r1 = rules.split_processed_bam.output.r1,
#r2 = rules.split_processed_bam.output.r2,
#
# sorted version 
#r1 = rules.sort_and_index_processed_bams.output.r1,
#r2 = rules.sort_and_index_processed_bams.output.r2,
#
# linedup version 
#r1 = rules.lineup_processed_bams.output.r1,
#r2 = rules.lineup_processed_bams.output.r2,
#
## Make perREfragStats for processed data
#rule make_perREfragStats_for_processed_alns_bioreps:
#    input:
#        r1 = rules.lineup_processed_bams.output.r1,
#        r2 = rules.lineup_processed_bams.output.r2,
#        frag = re_digestion_file,
#        fgc_map = re_fgc_map_file
#    output:
#        frag_stats = 'results/main/{cline}/hicnv/re_frag_stats/{cline}.{srr}.perREfragStats'
#    params:
#        merged_singletons = 'results/main/{cline}/hicnv/re_frag_stats/{cline}.{srr}.bwt2pairs.withSingles.mapq30.bam',
#        outdir = 'results/main/{cline}/hicnv/re_frag_stats/{cline}_allMap2FragmentsOutput',
#        prefix = '{cline}_allMap2FragmentsOutput'
#    resources:
#        nodes = 1,
#        ppn = 1,
#        mem_mb = 50000
#    log:
#        'results/main/{cline}/logs/rule_make_perREfragStats_for_4dn_reps_{cline}_{srr}.log'
#    shell:
#        r"""
#              # Merge singletons
#              # pysam will through some warnings about a missing index
#              # but this can be ignored according for forums
#              {config[python2]} workflow/scripts/mergeSAM-singletons.py \
#                            -f {input.r1} \
#                            -r {input.r2} \
#                            -o {params.merged_singletons} \
#                            -v \
#                            --single \
#                            -q 30 >> {log} 2>&1
#
#              mkdir -p {params.outdir}
#
#              # Map HiC fragments
#              {config[python2]} workflow/scripts/mapped_2hic_fragments.py \
#                            -f {input.frag} \
#                            -r {params.merged_singletons} \
#                            -s 100 \
#                            -l 800 \
#                            -d 1000 \
#                            --all \
#                            -v \
#                            -o {params.outdir} >> {log} 2>&1
#
#              # Sort the REfragStats
#              cat {params.outdir}/{wildcards.cline}.{wildcards.srr}.bwt2pairs.withSingles.mapq30.perREfragStats \
#                            | sort -k1,1 -k2,2n > {output} 2> {log}
#          """
#
## Filtering for chromosomes 1-22 + X.
## The main hicnv_v2.R script will not work because there are too
## few datapoints generated for chrM and chrY should also be excluded.
## Saved in the combined tech_combined folder to continue running with the 
## rest of the pipeline. After this rules combine_bio_reps_perREfragStats,
## and run_bio_hicnv_auto_bandwidth will complete running HiCNV.
#use rule filter_perREfragStats_for_main_chrs as filter_perREfragStats_for_processed_alns_bioreps with:
#    input:
#        frag_stats = rules.make_perREfragStats_for_processed_alns_bioreps.output.frag_stats,
#        chr_list = 'resources/chromosome_lists/main_chrs.txt'
#    output:
#        frag_stats = 'results/main/{cline}/hicnv/tech_combined/{cline}.{srr}.perREfragStats'
#    log:
#        'results/main/{cline}/logs/rule_filter_perREfragStats_for_4dn_reps_{cline}_{srr}.log'
#
#
#
#


# Run HiCnv after merging biological replicates
# I am also saving the output of HiCnv as a meta file
# because HiCnv is calculating the reference chromosome
rule run_bio_hicnv_auto_bandwidth_4dn: # (status: running)
    input:
        feat = re_fgc_map_sorted_chrflt_file,
        cov = rules.map_replicates_4dn.output
    params:
        parent_dir = 'results/main/{cline}/hicnv/bio_run_auto_bandwidth/',
        meta = '{cline}.meta.txt',
    output:
        directory('results/main/{cline}/hicnv/bio_run_auto_bandwidth/{cline}_hicnv')
    log:
        'results/main/{cline}/logs/rule_run_bio_hicnv_auto_bandwidth_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_run_bio_hicnv_auto_bandwidth_{cline}.bmk'
    resources:
        mem_mb = 50000,
        ppn = 1
    shadow: 'full'
    shell:
        r"""
            # extract the fragment cutoff from the config files
            re=$(grep {wildcards.cline} config/samplesheet.tsv | cut -f 2)
            fragcutoff=$(grep $re config/re_metadata.tsv | cut -f 3)
            echo "Restriction enzyme: ${{re}}" >> {log} 2>&1
            echo "Fragment cutoff: ${{fragcutoff}}" >> {log} 2>&1

            # run hicnv
            {config[R4]} workflow/scripts/hicnv_v3.R \
                --refeature={input.feat} \
                --coverage={input.cov} \
                --prefix={wildcards.cline} \
                --fragcutoff=$fragcutoff \
                --bandwidth=0 \
                --cpu {resources.ppn} | tee {params.meta} >> {log} 2>&1 \
                    || echo "Failed but saving output for debugging purposes." >> {log} 2>&1

            # making the outdir since snakemake doesn't premake it like it does for the
            # parent directory of an output file.
            mkdir -p {params.parent_dir} >> {log} 2>&1

            # moving the meta file into the hicnv results before the final move
            mkdir -p {wildcards.cline}_hicnv/logs/
            mv {params.meta} {wildcards.cline}_hicnv/logs/

            # move the HiCnv results to the run folder
            mv {wildcards.cline}_hicnv/ {output} >> {log} 2>&1
        """


ruleorder: run_bio_hicnv_auto_bandwidth > run_bio_hicnv_auto_bandwidth_4dn

# DEVELOPING
## Run HiCnv a chromosome at a time
##rules.make_perREfragStats.output.frag_stats
##rules.process_refeature.output.sorted_feat_map
## combined_frag_stats = 'results/main/{cline}/hicnv/combined/{cline}.{srr}.perREfragStats'
#rule run_tech_hicnv_auto_bandwidth_4dn:
#    input:
#        feat = re_fgc_map_sorted_chrflt_file,
#        cov = rules.combine_tech_reps_perREfragStats.output.combined_frag_stats,
#        bio_cns = 'results/main/{cline}/hicnv/bio_run_auto_bandwidth/{cline}_hicnv/{cline}.copyNumber.txt'
#    params:
#        parent_dir = 'results/main/{cline}/hicnv/tech_run_auto_bandwidth/',
#    output:
#        directory('results/main/{cline}/hicnv/tech_run_auto_bandwidth/{cline}_{srr}_hicnv')
#    log:
#        'results/main/{cline}/logs/rule_run_tech_hicnv_auto_bandwidth_{cline}_{srr}.log'
#    benchmark:
#        'results/main/{cline}/benchmarks/rule_run_tech_hicnv_auto_bandwidth_{cline}_{srr}.bmk'
#    resources:
#        mem_mb = 50000,
#        ppn = 1
#    shadow: 'full'
#    shell:
#        r"""
#            # extract the fragment cutoff from the config files
#            echo "extract the fragment cutoff from the config files" >> {log}
#            re=$(grep {wildcards.cline} config/samplesheet.tsv | cut -f 2)
#            fragcutoff=$(grep $re config/re_metadata.tsv | cut -f 3)
#            echo "Restriction enzyme: ${{re}}" >> {log} 2>&1
#            echo "Fragment cutoff: ${{fragcutoff}}" >> {log} 2>&1
#
#            # extract the reference chromosome from the corresponding biological rep
#            refchrom=$(tail -n 1 {input.bio_cns} | cut -f 8)
#
#            # run hicnv
#            echo "running hicnv" >> {log}
#            {config[R4]} workflow/scripts/hicnv_v3.R \
#                --refeature={input.feat} \
#                --coverage={input.cov} \
#                --prefix={wildcards.cline}.{wildcards.srr} \
#                --refchrom=$refchrom \
#                --fragcutoff=$fragcutoff \
#                --bandwidth=0 \
#                --cpu {resources.ppn} >> {log} 2>&1 || echo "Failed but saving output for debugging purposes." >> {log} 2>&1
#
#            # removing the pre-made snakemake outdir
#            echo "removing the pre-made snakemake outdir" >> {log}
#            mkdir -p {params.parent_dir} >> {log} 2>&1
#
#            # move the HiCnv results to the run folder
#            echo "move the HiCnv results to the run folder" >> {log}
#            mv {wildcards.cline}.{wildcards.srr}_hicnv/ {output} >> {log} 2>&1
#        """
