# Align the HiC data
# https://github.com/nservant/HiC-Pro
#rule oned_read_coverage:
#    input:
#        bam1 = rules.hicpro_align.output.bam1,
#        bam2 = rules.hicpro_align.output.bam2
#    output:
#        aln = 'data/{cline}/aln/{cline}.{srr}.sam',
#    log:
#        'logs/rule_hicpro_align_{cline}_{srr}.log'
#    shell:
#        """
#            scripts/Read_coverage_generation/run_1DReadCoverage.pl {input} {output} >> {log} 2>&1
#        """


# Download the mappability file for hg38 and do processing for HiCnv
rule download_hg38_mappability:
    params:
        map = 'refs/hg38_mappability/k50.Umap.MultiTrackMappability.bw',
        bedgraph = 'refs/hg38_mappability/k50.Umap.MultiTrackMappability.bedGraph'
    output:
        sorted_bg = 'refs/hg38_mappability/k50.Umap.MultiTrackMappability.sorted.bedGraph'
    shell:
        """
            ## Check the https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/ to download genome specific mappability files
            wget https://www.pmgenomics.ca/hoffmanlab/proj/bismap/trackhub/hg38/k50.Umap.MultiTrackMappability.bw \\
                -O {params.map}

            ## Process the bigWig file and create bedGraph file
            bigWigToBedGraph {params.map} {params.bedgraph}
            sort -k 1,1 -k2,2n {params.bedgraph} > {output}

            rm {params.bedgraph}
        """


# Process the feature file as specified in HiCnv
rule process_refeature:
    input:
        ref = 'refs/hg38/hg38.fa',
        dig = 'refs/restriction_enzymes/hg38_mboi_digestion.bed',
        map = 'refs/hg38_mappability/k50.Umap.MultiTrackMappability.sorted.bedGraph'
    params:
        extended = 'refs/restriction_enzymes/hg38_mboi_digestion.extended.bed',
        gc = 'refs/restriction_enzymes/hg38_mboi_digestion.extended.gc.bed',
        gc_map = 'refs/restriction_enzymes/hg38_mboi_digestion.extended.gc.map.bed',
        f_gc_map = 'refs/restriction_enzymes/hg38_mboi_digestion.extended.fragment.gc.map.bed'
    output:
        sorted_feat_map = 'refs/restriction_enzymes/hg38_mboi_digestion.extended.fragment.gc.map.sorted.bed'
    shell:
        r"""
            # Create extended 500bp restriction fragment file
            awk '{{print $1"\t"$2"\t"$2+250"\t"$4"\t"$2"\t"$3"\t"$3-$2"\n"$1"\t"$3-250"\t"$3"\t"$4"\t"$2"\t"$3"\t"$3-$2}}' {input.dig} \
                | awk '{{if($2 >= 0){{print}}}}'| sortBed > {params.extended}

            # Find GC percentage of 500bp regions
            bedtools nuc -fi {input.ref} -bed {params.extended} > {params.gc}

            # Map the mappability over GC content file
            bedtools map -a {params.gc} -b {input.map} -c 4 -o mean > {params.gc_map}

            # Create F_GC_MAP file
            perl scripts/F_GC_MAP_Files/gc_map_per_fragment.pl {params.gc_map} {input.dig} > {params.f_gc_map}

            # Sort the final F_GC_MAP file
            bedtools sort -i {params.f_gc_map} > {output}

            rm {params}
        """


# fgc_map file: HindIII_hg38.500.50.F_GC_MAP.bed
# Create this file as per your restriction fragment.
# Scripts to create this file are under ../scripts/F_GC_MAP_Files/ directory.
rule make_perREfragStats:
    input:
        fwd_bwt2_folder = 'data/A549/hicpro/bowtie_results/bwt2/data/$replicate_prefix\_1_*.bwt2merged.bam',
        rev_bwt2_folder = 'data/A549/hicpro/bowtie_results/bwt2/data/$replicate_prefix\_2_*.bwt2merged.bam',
        frag = rules.digest_reference_genome.output.mboi,
        fgc_map = rules.process_refeature.output.sorted_feat_map
    output:
        frag_stats = '' # {params.prefix}.perREfragStats
    params:
        rep_prefix = 'data/A549/hicpro/bowtie_results/bwt2/data/$replicate_prefix\_Rep_{i}',
        fwd_bwt2_folder = '',
        rev_bwt2_folder = ''
    resources:
        nodes = 1
        ppn = 2
        mem_mb = 50000
    log:
        'logs/rule_make_perREfragStats_{}.log'
    shell:
        """
              # Merge singletons
              python {config[hicnv_scripts]}/mergeSAM-singletons.py \
                            -f {params.fwd_bwt2_folder} \
                            -r {params.rwd_bwt2_folder} \
                            -o {params.prefix}.bwt2pairs.withSingles.mapq30.bam \
                            -v \
                            --single \
                            -q 30

              mkdir {params.prefix}\_allMap2FragmentsOutput

              # Map HiC fragments
              python {config[hicnv_scripts]}/mapped_2hic_fragments.py \
                            -f {input.frag} \
                            -r {params.prefix}.bwt2pairs.withSingles.mapq30.bam \
                            -s 100 \
                            -l 800 \
                            -d 1000 \
                            --all \
                            -v \
                            -o {params.prefix}\_allMap2FragmentsOutput

              # Sort the REfragStats
              cat {params.prefix}\_allMap2FragmentsOutput/{params.prefix}.bwt2pairs.withSingles.mapq30.perREfragStats \
                            | sort -k1,1 -k2,2n > {params.prefix}.perREfragStats
          """


## Run HiCnv a chromosome at a time
#rule run_hicnv:
#    input:
#        feat = rules.process_refeature.output.sorted_feat_map,
#        cov = rules.one_dim_3div_bedpe_hicnv_version.output
#    output:
#        'output/{dataset}/{cline}/copy_numbers/{chrom}.test'
#    log:
#        'logs/run_hicnv_{dataset}_{cline}_{chrom}.log'
#    shell:
#        """
#            {config[R4]} {config[hicnv]} \
#                --refeature={input.feat} \
#                --coverage={input.cov} \
#                --refchrom={wildcards.chrom} \
#                --prefix={wildcards.cline} >> {log} 2>&1
#        """
