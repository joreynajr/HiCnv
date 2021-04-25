cline="tixon"
fns="results/main/${cline}/reads/fastq_pair/ "
snakemake --profile workflow/profiles/local $@ $fns
