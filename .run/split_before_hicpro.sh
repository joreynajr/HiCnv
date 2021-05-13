fns=""
for cline in HeLa_MboI_G1_Snyc HeLa_Kyoto_MboI_G1_Snyc_STAG1_Depleted HeLa_Kyoto_MboI_G1_Snyc_STAG2_Depleted;
do
    fn="results/main/${cline}/reads/split_fastqs/"
    fns+="$fn "
done

echo "snakemake --profile workflow/profiles/pbs-torque/ $@ $fns"
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns

