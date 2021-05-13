# init the file names for snakemake output
fns=""

# NOT run under complete parallelization
clines="A549 HepG2 "
clines+="HL60S4_14_Pore HL60S4_5_Pore HL60S4"

# run under complete parallelization
clines="22Rv1 HeLa_Kyoto_MboI_G1_Snyc_STAG1_Depleted HeLa_Kyoto_MboI_G1_Snyc_STAG2_Depleted "
clines+="HeLa_MboI_G1_Snyc Nalm6"
#for cline in 22Rv1;
for cline in $clines;
do
    # print cline
    #echo "	cline: $cline"

    # get bowtie2 dir
    bwt2_dir="results/main/${cline}/hicpro_with_parallel/bowtie_results/bwt2/"

    for srr_dir in $(ls -d $bwt2_dir/*);
    do
        # get srr
        srr="$(basename $srr_dir)"

        # print srr id
        #echo "		srr: $srr"

        fn="results/main/${cline}/hicnv/combined/${cline}.${srr}.perREfragStats"
        fns+="$fn "
        break
    done
done
echo $fns

#snakemake --profile workflow/profiles/local $@ $fns
snakemake --profile workflow/profiles/pbs-torque $@ $fns
