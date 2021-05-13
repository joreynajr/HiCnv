cline="22Rv1"
srr="SRR7760384"
chr="chr1"
fn="results/main/${cline}/hicnv/run/${cline}_${srr}_hicnv/figures/${cline}.${srr}.${chr}.pdf"
snakemake --profile workflow/profiles/local $@ $fn
