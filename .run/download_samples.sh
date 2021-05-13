#for cline in A549 HepG2 22Rv1 HeLa_Kyoto_MboI_G1_Snyc_STAG1_Depleted HeLa_Kyoto_MboI_G1_Snyc_STAG2_Depleted HeLa_MboI_G1_Snyc HL60S4_14_Pore HL60S4_5_Pore HL60S4 Nalm6;
fns=""
for cline in A549 HepG2 HeLa_Kyoto_MboI_G1_Snyc_STAG1_Depleted HeLa_Kyoto_MboI_G1_Snyc_STAG2_Depleted HeLa_MboI_G1_Snyc HL60S4_14_Pore HL60S4_5_Pore;
do
    srrs=$(cat results/main/${cline}/reads/*.SRR_Acc_List.txt)
    for srr in $srrs;
    do
        # running the download only
        fn="results/main/${cline}/reads/${srr}_1.fastq.gz"

        # concat current file
        fns+="${fn} "
    done
done
echo "snakemake --profile workflow/profiles/pbs-torque/ $@ $fns"
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
