outputs=""
#for fn in $(ls data/*/sra/*.SRR_Acc_List.txt);
for fn in $(ls data/test/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    output="data/${cline}/hic_results/data/${srr}/${srr}.allValidPairs"
    outputs+="$output "
    break
done
echo $outputs
snakemake --profile  profiles/local $@ $outputs
#snakemake --profile profiles/pbs-torque $@ $outputs
