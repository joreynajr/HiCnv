outputs=""
#for fn in $(ls results/main/*/sra/*.SRR_Acc_List.txt);
for fn in $(ls results/main/test/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    output="results/main/${cline}/hicpro/hic_results/data/${srr}/${srr}_hg38.bwt2pairs.validPairs"
    outputs+="$output "
    break
done
echo $outputs
snakemake --profile  profiles/local $@ $outputs
#snakemake --profile profiles/pbs-torque $@ $outputs
