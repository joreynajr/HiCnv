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
