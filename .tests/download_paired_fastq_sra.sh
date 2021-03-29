outputs=""
#for fn in $(ls results/main/*/sra/*.SRR_Acc_List.txt);
for fn in $(ls results/main/A549/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    output="results/main/${cline}/sra/${srr}_1.fastq.gz"
    outputs="$output "
done
#snakemake --profile workflow/profiles/local $@ $outputs
snakemake --profile workflow/profiles/pbs-torque $@ $outputs
