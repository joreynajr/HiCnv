#snakemake --profile workflow/profiles/pbs-torque $@ results/refs/restriction_enzymes/hg38_mboi_digestion.bed
snakemake --profile workflow/profiles/local $@ results/refs/restriction_enzymes/hg38_mboi_digestion.bed
