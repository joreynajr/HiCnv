outputs=""
#for fn in $(ls data/*/sra/*.SRR_Acc_List.txt);
for fn in $(ls data/22Rv1/sra/*.SRR_Acc_List.txt);
do
    cline=$(basename $fn | cut -f 1 -d "." )
    srr=$(cat $fn)
    #output="data/${cline}/hicpro/${cline}.${srr}.1.bwt2merged.bam"
    output="data/${cline}/${srr}/hicpro/${cline}.${srr}.ran.flag"
    outputs+="$output "
done
#echo $outputs
#snakemake --profile $@ profiles/local $outputs
snakemake --profile $@ profiles/pbs-torque $outputs
