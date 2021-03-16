# rule to download sra files which derived from paired sequencing
# using the grabseqs tool and modified code from the hisss snakemake pipeline
# https://github.com/louiejtaylor/hisss/blob/master/rules/sra_paired.rules
rule download_paired_fastq_sra:
	output:
		r1 = 'data/sra/{cline]/{srr}_1.fastq.gz',
		r2 = 'data/sra/{cline]/[srr}_2.fastq.gz'
	params:
		outdir = '',
		sample = '{sample}'
	threads: 4
	shell:
		"""
		grabseqs sra -t {threads} -f -r 4 --no_parsing -o {params.outdir} {params.samp}
		"""
