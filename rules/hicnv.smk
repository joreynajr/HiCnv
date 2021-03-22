
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
        dig = 'refs/restriction_enzymes/hg38_{re}_digestion.bed',
        map = 'refs/hg38_mappability/k50.Umap.MultiTrackMappability.sorted.bedGraph'
    params:
        extended = 'refs/restriction_enzymes/hg38_{re}_digestion.extended.bed',
        gc = 'refs/restriction_enzymes/hg38_{re}_digestion.extended.gc.bed',
        gc_map = 'refs/restriction_enzymes/hg38_{re}_digestion.extended.gc.map.bed',
        f_gc_map = 'refs/restriction_enzymes/hg38_{re}_digestion.extended.fragment.gc.map.bed'
    output:
        sorted_feat_map = 'refs/restriction_enzymes/hg38_{re}_digestion.extended.fragment.gc.map.sorted.bed'
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


# Create this file as per your restriction fragment.
# Scripts to create this file are under ../scripts/F_GC_MAP_Files/ directory.
rule make_perREfragStats:
    input:
        fwd_alns = 'data/{cline}/hicpro/bowtie_results/bwt2/{srr}/{srr}_1_hg38.bwt2merged.bam',
        rev_alns = 'data/{cline}/hicpro/bowtie_results/bwt2/{srr}/{srr}_2_hg38.bwt2merged.bam',
        frag = 'refs/restriction_enzymes/hg38_hindiii_digestion.bed',
        fgc_map = 'refs/restriction_enzymes/hg38_hindiii_digestion.extended.fragment.gc.map.sorted.bed'
    output:
        frag_stats = 'data/{cline}/hicnv/{srr}.perREfragStats'
    params:
        merged_singletons = 'data/{cline}/hicnv/{srr}.bwt2pairs.withSingles.mapq30.bam',
        outdir = 'data/{cline}/hicnv/{srr}_allMap2FragmentsOutput',
        prefix = '{srr}_allMap2FragmentsOutput'
    resources:
        nodes = 1,
        ppn = 2,
        mem_mb = 50000
    log:
        'logs/rule_make_perREfragStats_{cline}_{srr}.log'
    shell:
        """
              # Merge singletons
              # pysam will through some warnings about a missing index
              # but this can be ignored according for forums
              {config[python2]} {config[hicnv_scripts]}/mergeSAM-singletons.py \
                            -f {input.fwd_alns} \
                            -r {input.rev_alns} \
                            -o {params.merged_singletons} \
                            -v \
                            --single \
                            -q 30

              mkdir {params.outdir}

              # Map HiC fragments
              {config[python2]} {config[hicnv_scripts]}/mapped_2hic_fragments.py \
                            -f {input.frag} \
                            -r {params.merged_singletons} \
                            -s 100 \
                            -l 800 \
                            -d 1000 \
                            --all \
                            -v \
                            -o {params.outdir}

              # Sort the REfragStats
              cat {params.outdir}/{wildcards.srr}.bwt2pairs.withSingles.mapq30.perREfragStats \
                            | sort -k1,1 -k2,2n > {output}
          """


# One dimensionalzie the read coverage data
rule onedim_read_coverage:
    input:
        bam1 = rules.hicpro_align_only.output.bam1,
        bam2 = rules.hicpro_align_only.output.bam2
    output:
        aln = 'data/{cline}/hicnv/{cline}.{srr}.onedim.bed',
    log:
        'logs/rule_onedim_read_coverage_{cline}_{srr}.log'
    shell:
        """
            perl scripts/Read_coverage_generation/run_1DReadCoverage.pl {input} {output} >> {log} 2>&1
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
