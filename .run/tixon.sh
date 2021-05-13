# mapping
cline="tixon"
srr="SRR4002640"

# align samples
#fn="results/main/${cline}/hicpro/bowtie_results/bwt2/${cline}/${srr}_1_hg38.bwt2merged.bam"

## merge samples
#fn="results/main/tixon/hicpro/hic_results/data/tixon/tixon.allValidPairs"

## build matrix
#fn="results/main/tixon/hicpro/hic_results/matrix/tixon/raw/"

## ice matrix
fn="results/main/${cline}/hicpro/hic_results/matrix/${cline}/iced/"

# running the whole pipeline
#fn="results/main/${cline}/hicnv/${cline}.${srr}_hicnv/CNV_Estimation/${cline}.${srr}.cnv.bedGraph"
#snakemake --profile workflow/profiles/pbs-torque/ $@ $fn
snakemake --profile workflow/profiles/local/ $@ $fn
