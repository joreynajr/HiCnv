outputs=""
#for fn in $(ls results/main/*/sra/*.SRR_Acc_List.txt);
for fn in $(ls results/main/test_copy/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    output="results/main/${cline}/hicpro/hic_results/matrix/${srr}/iced/10000/${srr}_10000_iced.matrix"
    outputs+="$output "
done
echo $outputs
snakemake --profile workflow/profiles/local $@ $outputs
#snakemake --profile workflow/profiles/pbs-torque $@ $outputs
