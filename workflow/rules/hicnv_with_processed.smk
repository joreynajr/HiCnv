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
