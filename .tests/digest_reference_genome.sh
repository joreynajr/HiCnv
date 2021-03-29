#snakemake --profile profiles/pbs-torque $@ refs/restriction_enzymes/hg38_mboi_digestion.bed
snakemake --profile profiles/local $@ refs/restriction_enzymes/hg38_mboi_digestion.bed
