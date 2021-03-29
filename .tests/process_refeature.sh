#snakemake --profile profiles/pbs-torque $@ refs/restriction_enzymes/hg38_mboi_digestion.extended.fragment.gc.map.sorted.bed
snakemake --profile profiles/local $@ refs/restriction_enzymes/hg38_mboi_digestion.extended.fragment.gc.map.sorted.bed
