outputs=""

cline="tixon"
output="results/main/${cline}/hicpro_with_adj/bowtie_results/bwt2/${cline}/"
outputs+="$output "

echo "snakemake --profile  workflow/profiles/local $@ $outputs"
snakemake --profile  workflow/profiles/local $@ $outputs
#snakemake --profile workflow/profiles/pbs-torque $@ $outputs
