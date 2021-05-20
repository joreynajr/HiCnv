# Run HiCnv after merging biological replicates
# I am also saving the output of HiCnv as a meta file
# because HiCnv is calculating the reference chromosome
rule run_bio_hicnv_auto_bandwidth:
    input:
        feat = re_fgc_map_sorted_chrflt_file,
        cov = rules.combine_bio_reps_perREfragStats.output.combined_frag_stats
    params:
        parent_dir = 'results/main/{cline}/hicnv/bio_run_auto_bandwidth/',
        meta = '{cline}.meta.txt',
        new_cov_format = 'results/main/{cline}/hicnv/bio_combined/{cline}.new_format.perREfragStats'
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


# Run HiCnv a chromosome at a time
#rules.make_perREfragStats.output.frag_stats
#rules.process_refeature.output.sorted_feat_map
# combined_frag_stats = 'results/main/{cline}/hicnv/combined/{cline}.{srr}.perREfragStats'
rule run_tech_hicnv_auto_bandwidth:
    input:
        feat = re_fgc_map_sorted_chrflt_file,
        cov = rules.combine_tech_reps_perREfragStats.output.combined_frag_stats,
        bio_cns = 'results/main/{cline}/hicnv/bio_run_auto_bandwidth/{cline}_hicnv/{cline}.copyNumber.txt'
    params:
        parent_dir = 'results/main/{cline}/hicnv/tech_run_auto_bandwidth/',
    output:
        directory('results/main/{cline}/hicnv/tech_run_auto_bandwidth/{cline}_{srr}_hicnv')
    log:
        'results/main/{cline}/logs/rule_run_tech_hicnv_auto_bandwidth_{cline}_{srr}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_run_tech_hicnv_auto_bandwidth_{cline}_{srr}.bmk'
    resources:
        mem_mb = 50000,
        ppn = 1
    shadow: 'full'
    shell:
        r"""
            # extract the fragment cutoff from the config files
            echo "extract the fragment cutoff from the config files" >> {log}
            re=$(grep {wildcards.cline} config/samplesheet.tsv | cut -f 2)
            fragcutoff=$(grep $re config/re_metadata.tsv | cut -f 3)
            echo "Restriction enzyme: ${{re}}" >> {log} 2>&1
            echo "Fragment cutoff: ${{fragcutoff}}" >> {log} 2>&1

            # extract the reference chromosome from the corresponding biological rep
            refchrom=$(tail -n 1 {input.bio_cns} | cut -f 8)

            # run hicnv
            echo "running hicnv" >> {log}
            {config[R4]} workflow/scripts/hicnv_v3.R \
                --refeature={input.feat} \
                --coverage={input.cov} \
                --prefix={wildcards.cline}.{wildcards.srr} \
                --refchrom=$refchrom \
                --fragcutoff=$fragcutoff \
                --bandwidth=0 \
                --cpu {resources.ppn} >> {log} 2>&1 || echo "Failed but saving output for debugging purposes." >> {log} 2>&1

            # removing the pre-made snakemake outdir
            echo "removing the pre-made snakemake outdir" >> {log}
            mkdir -p {params.parent_dir} >> {log} 2>&1

            # move the HiCnv results to the run folder
            echo "move the HiCnv results to the run folder" >> {log}
            mv {wildcards.cline}.{wildcards.srr}_hicnv/ {output} >> {log} 2>&1
        """


# Helper function to obtain all perREfragStats for a technical replicate
def get_bio_rep_auto_bw_cnv_est(wildcards):
    fn = 'results/main/{cline}/hicnv/bio_run_auto_bandwidth/{cline}_hicnv/CNV_Estimation/{cline}.chr*.cnv.bedGraph'
    input_fns = glob.glob(fn.format(**wildcards))
    return(sorted(input_fns))


# merged the cnv beds
rule merge_bio_hicnv_cnv_estimation_bedGraphs:
    input:
        get_bio_rep_auto_bw_cnv_est
    output:
        merged = 'results/main/{cline}/hicnv/bio_run_auto_bandwidth/{cline}_hicnv/figures/{cline}.bedGraph'
    shell:
        """
            cat {input} > {output}
        """


# plot for auto bandwidth using custom script
use rule plot_tech_hicnv_bedpe as plot_bio_hicnv_with_bw_bedpe with:
    input:
        bedgraph = rules.merge_bio_hicnv_cnv_estimation_bedGraphs.output.merged
    output:
        image = 'results/main/{cline}/hicnv/bio_run_auto_bandwidth/{cline}_hicnv/figures/{cline}.png'
    log:
        'results/main/{cline}/logs/rule_plot_bio_hicnv_auto_bandwidth_bedpe_{cline}.log'


# Helper function to obtain all perREfragStats for a technical replicate
def get_tech_rep_auto_bw_cnv_est(wildcards):
    fn = 'results/main/{cline}/hicnv/tech_run_auto_bandwidth/{cline}_{srr}_hicnv/CNV_Estimation/{cline}.{srr}.chr*.cnv.bedGraph'
    input_fns = glob.glob(fn.format(**wildcards))
    return(sorted(input_fns))


# merged the cnv beds
rule merge_tech_hicnv_cnv_estimation_bedGraphs:
    input:
        get_tech_rep_auto_bw_cnv_est
    output:
        merged = 'results/main/{cline}/hicnv/tech_run_auto_bandwidth/{cline}_{srr}_hicnv/figures/{cline}.{srr}.bedGraph'
    shell:
        """
            cat {input} > {output}
        """


# plot for auto bandwidth using custom script
use rule plot_tech_hicnv_bedpe as plot_tech_hicnv_with_bw_bedpe with:
    input:
        bedgraph = rules.merge_tech_hicnv_cnv_estimation_bedGraphs.output.merged
    output:
        image = 'results/main/{cline}/hicnv/tech_run_auto_bandwidth/{cline}_{srr}_hicnv/figures/{cline}.{srr}.png'
    log:
        'results/main/{cline}/logs/rule_plot_tech_hicnv_auto_bandwidth_bedpe_{cline}_{srr}.log'


# Helper function to obtain all perREfragStats for a technical replicate
def get_all_tech_rep_auto_bw_cnv_est(wildcards):
    fn = 'results/main/{cline}/hicnv/tech_run_auto_bandwidth/{cline}_*_hicnv/figures/{cline}.*.bedGraph'
    input_fns = glob.glob(fn.format(**wildcards))
    return(sorted(input_fns))

# plot all tech reps on a single plot for auto bandwidth using custom script
# output will be saved within the merged hicnv output
rule plot_all_tech_hicnv_bedpe:
    input:
        bedgraphs = get_all_tech_rep_auto_bw_cnv_est
    output:
        image = 'results/main/{cline}/hicnv/bio_run_auto_bandwidth/{cline}_hicnv/figures/{cline}.all_tech_reps.png'
    log:
        'results/main/{cline}/logs/rule_plot_all_tech_hicnv_bedpe_{cline}.log'
    params:
        max_cn = 11
    resources:
        mem_mb = 8000
    shell:
        r"""
            python workflow/scripts/plot_hicnv_bedgraphs_for_tech_reps.py \
                        --bedgraphs {input} \
                        --outfn {output} \
                        --max-cn {params.max_cn} >> {log} 2>&1
        """





























