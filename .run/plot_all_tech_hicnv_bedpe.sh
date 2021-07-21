# init the file names for snakemake output
fns=""

# completely done
clines="22Rv1 HeLa_Kyoto_MboI_G1_Snyc_STAG1_Depleted HeLa_Kyoto_MboI_G1_Snyc_STAG2_Depleted "
clines+="HeLa_MboI_G1_Snyc Nalm6 "
clines+="A549 HepG2 "
#clines+="HL60S4_14_Pore HL60S4_5_Pore HL60S4"

for cline in $clines;
do
    # print cline
    #echo "	cline: $cline"
    fn="results/main/${cline}/hicnv/bio_run_auto_bandwidth/${cline}_hicnv/figures/${cline}.all_tech_reps.png"
    fns+="$fn "
done
echo $fns

snakemake --profile workflow/profiles/local $@ $fns
#snakemake --profile workflow/profiles/pbs-torque $@ $fns
