outputs=""
#for fn in $(ls results/main/*/sra/*.SRR_Acc_List.txt);
#for fn in $(ls results/main/test/sra/*.SRR_Acc_List.txt);
#for fn in $(ls results/main/test_copy/sra/*.SRR_Acc_List.txt);
for fn in $(ls results/main/22Rv1/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    output="results/main/${cline}/hicnv/${cline}_${srr}_hicnv_final.test"
    outputs+="$output "
    break
done
echo $outputs
snakemake --profile workflow/profiles/local $@ $outputs
#snakemake --profile profiles/pbs-torque $@ $outputs
