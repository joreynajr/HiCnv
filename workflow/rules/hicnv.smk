# Download the mappability file for hg38 and do processing for HiCnv
rule download_hg38_mappability: # Restruct test done
    params:
        map = 'results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.bw',
        bedgraph = 'results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.bedGraph'
    output:
        sorted_bg = 'results/refs/hg38_mappability/k50.Umap.MultiTrackMappability.sorted.bedGraph'
    log:
        'logs/rule_download_hg38_mappability.log'
    shell:
        r"""
            ## Check the https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/ to download genome specific mappability files
            wget https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/hg38/k50.Umap.MultiTrackMappability.bw \
                -O {params.map} >> {log} 2>&1

            ## Process the bigWig file and create bedGraph file
            workflow/scripts/bigWigToBedGraph {params.map} {params.bedgraph} >> {log} 2>&1
            sort -k 1,1 -k2,2n {params.bedgraph} > {output} >> {log} 2>&1
            rm {params.bedgraph} >> {log} 2>&1
        """


# Process the feature file as specified in HiCnv
rule process_refeature: # Restruct test done
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
            perl scripts/F_GC_MAP_Files/gc_map_per_fragment.pl {params.gc_map} {input.dig} > {params.f_gc_map} 2> {log}

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
        fwd_alns = 'results/main/{cline}/hicpro/bowtie_results/bwt2/{cline}/{cline}_1_hg38.bwt2merged.bam',
        rev_alns = 'results/main/{cline}/hicpro/bowtie_results/bwt2/{cline}/{cline}_2_hg38.bwt2merged.bam',
        frag = re_digestion_file,
        fgc_map = re_fgc_map_file
    output:
        frag_stats = 'results/main/{cline}/hicnv/{cline}.perREfragStats'
    params:
        merged_singletons = 'results/main/{cline}/hicnv/{cline}.bwt2pairs.withSingles.mapq30.bam',
        outdir = 'results/main/{cline}/hicnv/{cline}_allMap2FragmentsOutput',
        prefix = '{cline}_allMap2FragmentsOutput'
    resources:
        nodes = 1,
        ppn = 2,
        mem_mb = 50000
    log:
        'results/main/{cline}/logs/rule_make_perREfragStats_{cline}.log'
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
              cat {params.outdir}/{wildcards.cline}.bwt2pairs.withSingles.mapq30.perREfragStats \
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
        frag_stats = 'results/main/{cline}/hicnv/{cline}.chrflt.perREfragStats'
    log:
        'results/main/{cline}/logs/rule_filter_perREfragStats_for_main_chrs_{cline}.log'
    shell:
        r"""
            python workflow/scripts/filter_chrs_from_bed.py \
                --bed-like {input.frag_stats} \
                --include-list {input.chr_list} \
                -o {output} >> {log} 2>&1
        """


# get the restriction enzyme digestion files
def re_fgc_map_sorted_chrflt_file(wildcards):
    re = SAMPLESHEET.loc[wildcards.cline, 're']
    config = 'results/refs/restriction_enzymes/hg38_{}_digestion.extended.fragment.gc.map.sorted.chrflt.bed'.format(re)
    return(config)


# Run HiCnv a chromosome at a time
#rules.make_perREfragStats.output.frag_stats
#rules.process_refeature.output.sorted_feat_map
rule run_hicnv:
    input:
        feat = re_fgc_map_sorted_chrflt_file,
        cov = rules.filter_perREfragStats_for_main_chrs.output.frag_stats
    params:
        outdir = directory('results/main/{cline}/hicnv/')
    output:
        'results/main/{cline}/hicnv/{cline}.{cline}_hicnv/CNV_Estimation/{cline}.{cline}.cnv.bedGraph'
    log:
        'results/main/{cline}/logs/rule_run_hicnv_{cline}_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_run_hicnv_{cline}_{cline}.bmk'
    resources:
        ppn = 4
    #shadow: 'minimal'
    shell:
        r"""
            {config[R4]} workflow/scripts/hicnv_v2.R \
                --refeature={input.feat} \
                --coverage={input.cov} \
                --prefix={wildcards.cline} \
                --fragcutoff=150 \
                --refchrom=chr1 \
                --cpu {resources.ppn} >> {log} 2>&1

            mkdir -p {params.outdir} >> {log} 2>&1
            mv {wildcards.cline}_hicnv/ {params.outdir} >> {log} 2>&1
        """
