fns=""

# technical rep 1
cline="LNCaP"
srr="4DNESNHN919R-B1-T1"
fns+="results/main/${cline}/hicnv/re_frag_stats/${cline}.${srr}.perREfragStats "

# technical rep 2
cline="LNCaP"
srr="4DNESNHN919R-B1-T2"
#fns+="results/main/${cline}/hicnv/re_frag_stats/${cline}.${srr}.perREfragStats "

#snakemake --profile workflow/profiles/local $@ $fns
snakemake --profile workflow/profiles/pbs-torque/ $@ $fns
