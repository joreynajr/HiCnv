#snakemake --profile workflow/profiles/pbs-torque $@ results/refs/restriction_enzymes/hg38_mboi_digestion.extended.fragment.gc.map.sorted.bed
snakemake --profile workflow/profiles/local $@ results/refs/restriction_enzymes/hg38_mboi_digestion.extended.fragment.gc.map.sorted.bed
