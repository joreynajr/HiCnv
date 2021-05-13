cline='T47D'
libs=$(ls results/main/T47D/reads/T47D.ENC*.tsv | sed 's!.*/!!' | cut -f 2 -d ".")

fns=""
for lib in $libs;
do
    # running the download only 
    #fn="results/main/${cline}/sra/${lib}_1.fastq.gz"

    # mapping
    #fn="results/main/${cline}/hicpro/bowtie_results/bwt2/${lib}/${lib}_1_hg38.bwt2merged.bam"

    # running the whole pipeline
    fn="results/main/${cline}/hicnv/${cline}.${lib}_hicnv/CNV_Estimation/${cline}.${lib}.cnv.bedGraph"
    fns+=" ${fn}"
    break
done
echo "snakemake --profile workflow/profiles/pbs-torque/ $@ $fns"
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
