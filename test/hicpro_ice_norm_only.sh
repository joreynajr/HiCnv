outputs=""
#for fn in $(ls data/*/sra/*.SRR_Acc_List.txt);
for fn in $(ls data/test/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    output="data/${cline}/hicpro/hic_results/matrix/${srr}/iced/10000/${srr}_10000_iced.matrix"
    outputs+="$output "
    break
done
echo $outputs
snakemake --profile  profiles/local $@ $outputs
#snakemake --profile profiles/pbs-torque $@ $outputs
