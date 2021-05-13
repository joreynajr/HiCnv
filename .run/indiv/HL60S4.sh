cline='HL60S4'
srrs=$(cat results/main/${cline}/sra/*.SRR_Acc_List.txt)
fns=""
for srr in $srrs;
do
    # running the download only
    fn="results/main/${cline}/sra/${srr}_1.fastq.gz"

    # running the whole pipeline
    # fn="results/main/${cline}/hicnv/${cline}.${srr}_hicnv/CNV_Estimation/${cline}.${srr}.cnv.bedGraph"
    fns+=" ${fn}"
done
echo "snakemake --profile workflow/profiles/pbs-torque/ $@ $fns"
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
#snakemake --profile workflow/profiles/local/ $@ $fns
