cline="HL60S4_14_Pore"
fns="results/main/${cline}/hicpro/renamed_fastqs/${cline}/"
echo "snakemake --profile  workflow/profiles/pbs-torque $@ $fns"
snakemake --profile workflow/profiles/pbs-torque $@ $fns
