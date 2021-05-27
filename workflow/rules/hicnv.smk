#####################################################################
############## Methods developed NOT using auto bandwidth ###########
#####################################################################

# Download the mappability file for hg38 and do processing for HiCnv
rule download_hg38_mappability:
    params:
        map = 'results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.bw',
        bedgraph = 'results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.bedGraph'
    output:
        sorted_bg = 'results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.sorted.bedGraph'
    resources:
        mem_mb = 16000
    log:
        'results/refs/hg38_mappability/logs/rule_download_hg38_mappability.log'
    shell:
        r"""
            ## Check the https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/ to download genome specific mappability files
            wget https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/hg38/k50.Umap.MultiTrackMappability.bw \
                -O {params.map} >> {log} 2>&1

            ## Process the bigWig file and create bedGraph file
            workflow/scripts/bigWigToBedGraph {params.map} {params.bedgraph} >> {log} 2>&1
            sort -k 1,1 -k2,2n {params.bedgraph} > {output} 2>> {log}
            #rm {params.bedgraph} >> {log} 2>&1
        """


# Process the feature file as specified in HiCnv
rule process_refeature:
    input:
        ref = 'results/refs/hg38/hg38.fa',
        dig = 'results/refs/restriction_enzymes/hg38_{re}_digestion.bed',
        map = 'results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.sorted.bedGraph'
    params:
        extended = 'results/refs/restriction_enzymes/hg38_{re}_digestion.extended.bed',
        gc = 'results/refs/restriction_enzymes/hg38_{re}_digestion.extended.gc.bed',
        gc_map = 'results/refs/restriction_enzymes/hg38_{re}_digestion.extended.gc.map.bed',
        f_gc_map = 'results/refs/restriction_enzymes/hg38_{re}_digestion.extended.fragment.gc.map.bed'
    output:
        sorted_feat_map = 'results/refs/restriction_enzymes/hg38_{re}_digestion.extended.fragment.gc.map.sorted.bed'
    resources:
        mem_mb = 30000
    log:
        'results/refs/restriction_enzymes/logs/rule_process_refeature_{re}.log'
    shell:
        r"""
            # Create extended 500bp restriction fragment file
            awk '{{print $1"\t"$2"\t"$2+250"\t"$4"\t"$2"\t"$3"\t"$3-$2"\n"$1"\t"$3-250"\t"$3"\t"$4"\t"$2"\t"$3"\t"$3-$2}}' {input.dig} \
                | awk '{{if($2 >= 0){{print}}}}'| sortBed > {params.extended} 2> {log}

            # Find GC percentage of 500bp regions
            bedtools nuc -fi {input.ref} -bed {params.extended} > {params.gc} 2> {log}

            # Map the mappability over to the GC content file
            bedtools map -a {params.gc} -b {input.map} -c 4 -o mean > {params.gc_map} 2> {log}

            # Create F_GC_MAP file
            perl workflow/scripts/F_GC_MAP_Files/gc_map_per_fragment.pl {params.gc_map} {input.dig} > {params.f_gc_map} 2> {log}

            # Sort the final F_GC_MAP file
            bedtools sort -i {params.f_gc_map} > {output} 2> {log}

            rm {params}
        """


# Filtering for chromosomes 1-22 + X.
# The main hicnv_v2.R script will not work because there are too
# few datapoints generated for chrM and chrY should also be excluded.
rule filter_refeature_for_main_chrs:
    input:
        refeat = rules.process_refeature.output.sorted_feat_map,
        chr_list = 'resources/chromosome_lists/main_chrs.txt'
    output:
        refeat = 'results/refs/restriction_enzymes/hg38_{re}_digestion.extended.fragment.gc.map.sorted.chrflt.bed'
    log:
        'results/refs/restriction_enzymes/logs/rule_filter_refeature_for_main_chrs_{re}.log'
    shell:
        r"""
            python workflow/scripts/filter_chrs_from_bed.py \
                --bed-like {input.refeat} \
                --include-list {input.chr_list} \
                -o {output} >> {log} 2>&1
        """


# get the restriction enzyme digestion files
def re_fgc_map_file(wildcards):
    re = SAMPLESHEET.loc[wildcards.cline, 're']
    config = 'results/refs/restriction_enzymes/hg38_{}_digestion.extended.fragment.gc.map.sorted.bed'.format(re)
    return(config)


# Create this file as per your restriction fragment.
# Scripts to create this file are under ../scripts/F_GC_MAP_Files/ directory.
rule make_perREfragStats:
    input:
        fwd_alns = 'results/main/{cline}/hicpro_with_parallel/bowtie_results/bwt2/{srr}/{split}_{srr}_R1_hg38.bwt2merged.bam',
        rev_alns = 'results/main/{cline}/hicpro_with_parallel/bowtie_results/bwt2/{srr}/{split}_{srr}_R2_hg38.bwt2merged.bam',
        frag = re_digestion_file,
        fgc_map = re_fgc_map_file
    output:
        frag_stats = 'results/main/{cline}/hicnv/split_re_frag_stats/{cline}.{srr}.{split}.perREfragStats'
    params:
        merged_singletons = 'results/main/{cline}/hicnv/split_re_frag_stats/{cline}.{srr}.{split}.bwt2pairs.withSingles.mapq30.bam',
        outdir = 'results/main/{cline}/hicnv/split_re_frag_stats/{cline}_allMap2FragmentsOutput',
        prefix = '{cline}_allMap2FragmentsOutput'
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 50000
    log:
        'results/main/{cline}/logs/rule_make_perREfragStats_{cline}_{srr}_{split}.log'
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
              cat {params.outdir}/{wildcards.cline}.{wildcards.srr}.{wildcards.split}.bwt2pairs.withSingles.mapq30.perREfragStats \
                            | sort -k1,1 -k2,2n > {output} 2> {log}
          """


# Filtering for chromosomes 1-22 + X.
# The main hicnv_v2.R script will not work because there are too
# few datapoints generated for chrM and chrY should also be excluded.
rule filter_perREfragStats_for_main_chrs:
    input:
        frag_stats = rules.make_perREfragStats.output.frag_stats,
        chr_list = 'resources/chromosome_lists/main_chrs.txt'
    output:
        frag_stats = 'results/main/{cline}/hicnv/split_re_frag_stats/{cline}.{srr}.{split}.chrflt.perREfragStats'
    log:
        'results/main/{cline}/logs/rule_filter_perREfragStats_for_main_chrs_{cline}_{srr}_{split}.log'
    shell:
        r"""
            python workflow/scripts/filter_chrs_from_bed.py \
                --bed-like {input.frag_stats} \
                --include-list {input.chr_list} \
                -o {output} >> {log} 2>&1
        """


# Helper function to obtain all perREfragStats for a technical replicate
def get_tech_rep_chrflt_perREfragStats(wildcards):
    input_fns = []

    for fn in glob.glob('results/main/{cline}/hicpro_with_parallel/bowtie_results/bwt2/{srr}/*_{srr}_R1_hg38.bwt2merged.bam'.format(**wildcards)):
        split = os.path.basename(fn).split('_')[0]
        new_input = 'results/main/{cline}/hicnv/split_re_frag_stats/{cline}.{srr}.{split}.chrflt.perREfragStats'.format(split=split, **wildcards)
        input_fns.append(new_input)
    return(sorted(input_fns))


# Create this file as per your restriction fragment.
# Scripts to create this file are under ../scripts/F_GC_MAP_Files/ directory.
rule combine_tech_reps_perREfragStats:
    input:
        get_tech_rep_chrflt_perREfragStats
    params:
        re_fs_list = 'results/main/{cline}/hicnv/tech_combined/{cline}.{srr}.REfragStats.list.txt'
    output:
        combined_frag_stats = 'results/main/{cline}/hicnv/tech_combined/{cline}.{srr}.perREfragStats'
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 16000
    log:
        'results/main/{cline}/logs/rule_combine_tech_reps_perREfragStats_{cline}_{srr}.log'
    shell:
        r"""
            # separate the inputs by new line character
            for fn in {input};
            do
                echo $fn >> {params.re_fs_list}
            done

            # run the combine script
            perl workflow/scripts/combine_multiple_REfragStats.pl {params.re_fs_list} {output} >> {log} 2>&1
        """


# get the restriction enzyme digestion files
def re_fgc_map_sorted_chrflt_file(wildcards):
    re = SAMPLESHEET.loc[wildcards.cline, 're']
    config = 'results/refs/restriction_enzymes/hg38_{}_digestion.extended.fragment.gc.map.sorted.chrflt.bed'.format(re)
    return(config)


# Run HiCnv a chromosome at a time
#rules.make_perREfragStats.output.frag_stats
#rules.process_refeature.output.sorted_feat_map
# combined_frag_stats = 'results/main/{cline}/hicnv/combined/{cline}.{srr}.perREfragStats'
rule run_tech_hicnv:
    input:
        feat = re_fgc_map_sorted_chrflt_file,
        cov = rules.combine_tech_reps_perREfragStats.output.combined_frag_stats
    params:
        parent_dir = 'results/main/{cline}/hicnv/tech_run/',
    output:
        directory('results/main/{cline}/hicnv/tech_run/{cline}_{srr}_hicnv')
    log:
        'results/main/{cline}/logs/rule_run_tech_hicnv_{cline}_{srr}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_run_tech_hicnv_{cline}_{srr}.bmk'
    resources:
        mem_mb = 50000,
        ppn = 2
    shadow: 'full'
    shell:
        r"""
            # extract the fragment cutoff from the config files
            re=$(grep {wildcards.cline} config/samplesheet.tsv | cut -f 2)
            fragcutoff=$(grep $re config/re_metadata.tsv | cut -f 3)
            echo "Restriction enzyme: ${{re}}" >> {log} 2>&1
            echo "Fragment cutoff: ${{fragcutoff}}" >> {log} 2>&1

            # run hicnv
            {config[R4]} workflow/scripts/hicnv_v2.R \
                --refeature={input.feat} \
                --coverage={input.cov} \
                --prefix={wildcards.cline}.{wildcards.srr} \
                --fragcutoff=$fragcutoff \
                --cpu {resources.ppn} >> {log} 2>&1 || echo "Failed but saving output for debugging purposes." >> {log} 2>&1

            # removing the pre-made snakemake outdir
            mkdir -p {params.parent_dir} >> {log} 2>&1

            # move the HiCnv results to the run folder
            mv {wildcards.cline}.{wildcards.srr}_hicnv/ {output} >> {log} 2>&1
        """


# make the ini file required for plotting
# currently not put into use due to custom script
rule make_pyGenomeTracks_ini:
    input:
        bedgraph = 'results/main/{cline}/hicnv/run/{cline}_{srr}_hicnv/CNV_Estimation/{cline}.{srr}.cnv.bedGraph'
    output:
        ini = 'results/main/{cline}/hicnv/run/{cline}_{srr}_hicnv/figures/{cline}.{srr}.ini'
    log:
        'results/main/{cline}/logs/rule_make_pyGenomeTracks_{cline}_{srr}.log'
    shell:
        r"""
            {config[pyGenomeTracks]}/make_tracks_file --trackFiles {input.bedgraph} -o {output.ini}
        """


# plot using the ini file
# currently not put into use due to custom script
rule plot_with_pyGenomeTracks:
    input:
        ini = rules.make_pyGenomeTracks_ini.output.ini
    output:
        image = 'results/main/{cline}/hicnv/run/{cline}_{srr}_hicnv/figures/{cline}.{srr}.{chr}.pdf'
    log:
        'results/main/{cline}/logs/rule_plot_with_pyGenomeTracks_{cline}_{srr}_{chr}.log'
    shell:
        r"""
            {config[pyGenomeTracks]}/pyGenomeTracks --tracks {input.ini} --region {wildcards.chr}:0-10000000 --outFileName {output.image}
        """







# plot using custom script
rule plot_tech_hicnv_bedpe:
    input:
        bedgraph = 'results/main/{cline}/hicnv/tech_run/{cline}_{srr}_hicnv/CNV_Estimation/{cline}.{srr}.cnv.bedGraph'
    output:
        image = 'results/main/{cline}/hicnv/tech_run/{cline}_{srr}_hicnv/figures/{cline}.{srr}.coverage.png'
    log:
        'results/main/{cline}/logs/rule_plot_tech_hicnv_bedpe_{cline}_{srr}.log'
    params:
        max_cn = 11
    resources:
        mem_mb = 8000
    shell:
        r"""
            python workflow/scripts/plot_hicnv_bedgraph.py \
                        --bedgraph {input} \
                        --outfn {output} \
                        --max-cn {params.max_cn} >> {log} 2>&1
        """








# Helper function to obtain all perREfragStats for a technical replicate
def get_bio_rep_chrflt_perREfragStats(wildcards):
    input_fns = glob.glob('results/main/{cline}/hicnv/tech_combined/{cline}.*.perREfragStats'.format(**wildcards))
    return(sorted(input_fns))


# Create this file as per your restriction fragment.
# Scripts to create this file are under ../scripts/F_GC_MAP_Files/ directory.
rule combine_bio_reps_perREfragStats:
    input:
        get_bio_rep_chrflt_perREfragStats
    params:
        re_fs_list = 'results/main/{cline}/hicnv/bio_combined/{cline}.REfragStats.list.txt'
    output:
        combined_frag_stats = 'results/main/{cline}/hicnv/bio_combined/{cline}.perREfragStats'
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 16000
    log:
        'results/main/{cline}/logs/rule_combine_bio_reps_perREfragStats_{cline}.log'
    shell:
        r"""
            # separate the inputs by new line character
            num_files=$(echo {input} | wc -w)
            if [[ $num_files == 1 ]]; then
                ln {input} {output}
            else
                for fn in {input};
                do
                    echo $fn >> {params.re_fs_list}
                done

                # run the combine script
                perl workflow/scripts/combine_multiple_REfragStats.pl {params.re_fs_list} {output} >> {log} 2>&1
            fi
        """


# Run HiCnv after merging biological replicates
# I am also saving the output of HiCnv as a meta file
# because HiCnv is calculating the reference chromosome
rule run_bio_hicnv:
    input:
        feat = re_fgc_map_sorted_chrflt_file,
        cov = rules.combine_bio_reps_perREfragStats.output.combined_frag_stats
    params:
        parent_dir = 'results/main/{cline}/hicnv/bio_run/',
        meta = '{cline}.meta.txt'
    output:
        directory('results/main/{cline}/hicnv/bio_run/{cline}_hicnv')
    log:
        'results/main/{cline}/logs/rule_run_bio_hicnv_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_run_bio_hicnv_{cline}.bmk'
    resources:
        mem_mb = 50000,
        ppn = 2
    shadow: 'full'
    shell:
        r"""
            # extract the fragment cutoff from the config files
            re=$(grep {wildcards.cline} config/samplesheet.tsv | cut -f 2)
            fragcutoff=$(grep $re config/re_metadata.tsv | cut -f 3)
            echo "Restriction enzyme: ${{re}}" >> {log} 2>&1
            echo "Fragment cutoff: ${{fragcutoff}}" >> {log} 2>&1

            # run hicnv
            {config[R4]} workflow/scripts/hicnv_v2.R \
                --refeature={input.feat} \
                --coverage={input.cov} \
                --prefix={wildcards.cline} \
                --fragcutoff=$fragcutoff \
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

#####################################################################
############## Methods developed TO USE auto bandwidth ##############
#####################################################################

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
    log:
        'results/main/{cline}/logs/rule_merge_bio_hicnv_cnv_estimation_bedGraphs_{cline}.log'
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
    log:
        'results/main/{cline}/logs/rule_merge_tech_hicnv_cnv_estimation_bedGraphs_{cline}_{srr}.log'
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





























