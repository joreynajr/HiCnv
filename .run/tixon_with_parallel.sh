# mapping
cline="tixon"
srr="SRR4002640"

# align samples
#fn="results/main/${cline}/hicpro_with_parallel/bowtie_results/bowtie.running"
fn="results/main/${cline}/hicpro_with_parallel/bowtie_results/all.running"

# quality control
#fn="results/main/${cline}/hicpro_with_parallel/hic_results/pic/${cline}"

# processing
#fn="results/main/${cline}/hicpro_with_parallel/hic_results/data/${cline}/process.running"
#fn="results/main/${cline}/hicpro_with_parallel/hic_results/data/${cline}/process.running"

## merge samples
#fn="results/main/${cline}/hicpro_with_parallel/hic_results/data/${cline}/merge.running"

## build matrix
#fn="results/main/tixon/hicpro/hic_results/matrix/tixon/raw/"

## ice matrix
#fn="results/main/${cline}/hicpro/hic_results/matrix/${cline}/iced/"

# running the whole pipeline
#fn="results/main/${cline}/hicnv/${cline}.${srr}_hicnv/CNV_Estimation/${cline}.${srr}.cnv.bedGraph"
#snakemake --profile workflow/profiles/pbs-torque/ $@ $fn
snakemake --profile workflow/profiles/local/ $@ $fn
