# init the file names for snakemake output
fns=""

# completely done
clines="UCD52CR UCD52PR"
srr="custom"
for cline in $clines;
do
    # print cline
    #echo "	cline: $cline"

    # run HiCnv main
    #fn="results/main/${cline}/hicnv/run/${cline}_${srr}_hicnv/"

    # run HiCnv plotting
    fn="results/main/${cline}/hicnv/bio_run_auto_bandwidth/${cline}_hicnv"
    #fn="results/main/${cline}/hicnv/tech_run_auto_bandwidth/${cline}_${srr}_hicnv"
    fns+="$fn "
done
echo $fns

#cline="HL60S4_14_Pore"
#srr="SRR7291411"
#fns="results/main/${cline}/hicnv/run/${cline}_${srr}_hicnv/"

#snakemake --profile workflow/profiles/local $@ $fns
snakemake --profile workflow/profiles/pbs-torque $@ $fns
