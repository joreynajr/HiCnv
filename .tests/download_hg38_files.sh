#snakemake --profile workflow/profiles/pbs-torque $@ results/refs/hg38/hg38.fa.gz
snakemake --profile workflow/profiles/local $@ results/refs/hg38/hg38.fa.gz
