cline='Nalm6'
srrs=$(cat results/main/${cline}/reads/*.SRR_Acc_List.txt)
fns=""
for srr in $srrs;
do
    # running the download only 
    #fn="results/main/${cline}/sra/${srr}_1.fastq.gz"

    # mapping
    fn="results/main/${cline}/hicpro/bowtie_results/bwt2/${srr}/${srr}_1_hg38.bwt2merged.bam"

    # running the whole pipeline
    #fn="results/main/${cline}/hicnv/${cline}.${srr}_hicnv/CNV_Estimation/${cline}.${srr}.cnv.bedGraph"
    fns+=" ${fn}"
done
echo "snakemake --profile workflow/profiles/pbs-torque/ $@ $fns"
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
