fns=""

# technical rep 1
cline="LNCaP"
srr="4DNESNHN919R-B1-T1"
fns+="results/main/${cline}/4d_nucleome/${srr}/${srr}_R1.sorted.linedup.bam "

# technical rep 2
cline="LNCaP"
srr="4DNESNHN919R-B1-T2"
#fns+="results/main/${cline}/4d_nucleome/${srr}_R1.sorted.linedup.bam "

#snakemake --profile workflow/profiles/local $@ $fns
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
