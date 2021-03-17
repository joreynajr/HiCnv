outputs=""
#for fn in $(ls data/*/sra/*.SRR_Acc_List.txt);
for fn in $(ls data/A549/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    output="data/${cline}/sra/${srr}_1.fastq.gz"
    outputs="$output "
done
#snakemake --profile $@ profiles/local $outputs
snakemake --profile $@ profiles/pbs-torque $outputs
