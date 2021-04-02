outputs=""
for fn in $(ls results/main/test_copy/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    output="results/main/${cline}/hicnv/${cline}.${srr}_hicnv/CNV_Estimation/${cline}.${srr}.cnv.bedGraph"
    outputs+="$output "
    break
done

outputs2=""
##for fn in $(ls results/main/*/sra/*.SRR_Acc_List.txt);
#for fn in $(ls results/main/test_copy/sra/*.SRR_Acc_List.txt);
#do
#    cline=$(basename $fn | cut -f 1 -d "." )
#    srr=$(cat $fn)
#    output="results/main/${cline}/hicpro/hic_results/matrix/${srr}/iced/10000/${srr}_10000_iced.matrix"
#    outputs2+="$output "
#    break
#done

snakemake --profile workflow/profiles/local -n -f --dag $outputs $outputs2 | dot -Tsvg > images/graph.svg
