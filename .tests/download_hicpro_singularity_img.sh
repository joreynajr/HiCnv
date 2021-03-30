#snakemake --profile workflow/profiles/pbs-torque $@ resources/software/hicpro_latest_ubuntu.img
snakemake --profile workflow/profiles/local $@ resources/software/hicpro_latest_ubuntu.img
