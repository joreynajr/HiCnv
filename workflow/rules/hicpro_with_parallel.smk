# Helper function for to rule split_before_hicpro to obtain all SRR's 
# for a given cell line.
# Make sure to download the accession list from SRA
def get_r1_r2_fastqs_for_splitting(wildcards):

    # init lists to collect r1 and r2 files
    r1s = []
    r2s = []

    # load the fastq meta data 
    fastq_sra = read_table('config/fastq_sra_meta.tsv')
    fastq_4dn = read_table('config/fastq_4dn_meta.tsv')

    # checking for the current cell line in both list
    # each cell line should only come from a single list OR 
    # the cell line can be given a slightly different name.
    if (wildcards.cline in fastq_sra.cline.tolist()) and (wildcards.cline in fastq_4dn.cline.tolist()):
        msg = 'Cell line is in both config/fastq_sra_meta.tsv and '
        msg += 'config/fastq_sra_meta.tsv. Please remove cell line '
        msg += 'from one fastq config file before proceeding or give '
        msg += 'a different name. Pipeline terminaled without completion.'
        raise Exception(msg)

    # handle the r1 and r2 files when sra is the source
    if wildcards.cline in fastq_sra.cline.tolist():
        final_sra = fastq_sra[(fastq_sra.cline == wildcards.cline) & (fastq_sra.file_type == 'fastq')]
        for acc in sorted(final_sra.exp_acc.tolist()):
            r1 = 'results/main/{cline}/reads/{acc}_1.fastq.gz'.format(cline=wildcards.cline, acc=acc)
            r1s.append(r1)
            r2 = 'results/main/{cline}/reads/{acc}_2.fastq.gz'.format(cline=wildcards.cline, acc=acc)
            r2s.append(r2)
            
    # handle the r1 and r2 files when 4DN is the source
    elif wildcards.cline in fastq_4dn.cline.tolist():
        # filter for data from this cell line
        final_4dn = fastq_4dn[(fastq_4dn.cline == wildcards.cline) & (fastq_4dn.file_type == 'fastq')]

        # get the r1 and r2 names
        for i, sr in final_4dn.iterrows():
    
            # created the accession id
            acc = '{}-B{}-T{}'.format(sr.exp_acc, sr.biorep, sr.techrep)

            # name the fastq
            fastq = 'results/main/{cline}/reads/{acc}_{rnum}.fastq.gz'
            fastq = fastq.format(cline=wildcards.cline, acc=acc, rnum=sr.r1_or_r2)

            # add r1 or r2 to the corresponding list
            if sr.r1_or_r2 in [1, '1']:
                r1s.append(fastq)
            elif sr.r1_or_r2 in [2, '2']:
                r2s.append(fastq)
            else:
                msg = 'cline: {}, exp_acc: {} has an incorrect r1_or_r2: {}. '
                msg += 'Pipeline terminaled without completion.'
                msg = msg.format(wildcards.cline, acc, sr.r1_or_r2)
                raise Exception(msg)

    # sample is in neither fastq samplesheet
    else:
        msg = 'Cell line is not in config/fastq_sra_meta.tsv nor '
        msg += 'config/fastq_4dn_meta.tsv. Please add before proceeding. '
        msg += 'Pipeline terminaled without completion.'
        raise Exception(msg)

    d = {'r1s': r1s, 'r2s': r2s}
    return(d)


# Splitting the original R1 and R2 fastq's because they are really big
# Current this is set up to run serially but by adding an srr wildcard
# I can set it up to run in a more parallel way. The idea right now
# is that you would have a single large >70gb file per sample (bio-rep).
rule split_before_hicpro:
    input:
        unpack(get_r1_r2_fastqs_for_splitting)
    output:
        outdir = directory('results/main/{cline}/reads/split_fastqs/'),
        split_complete = touch('results/main/{cline}/reads/split_fastqs/split.complete')
    params:
        nreads = 50000000
    log:
        'results/main/{cline}/logs/rule_split_before_hicpro_{cline}.log'
    resources:
        nodes = 1,
        ppn = 1,
        mem_mb = 8000
    shell:
        r"""
            mkdir -p {output.outdir}

            # splitting the R1's
            for fq in {input.r1s} {input.r2s};
            do
                {config[python_hicpro]} {config[hicpro_dir]}/bin/utils/split_reads.py \
                                    --results_folder {output.outdir} \
                                    --nreads {params.nreads} \
                                    ${{fq}}
            done
        """


# Helper function to obtain all for a given cell line
# Make sure to download the accession list from SRA
def get_split_r1_r2_fastqs(wildcards):

    # init lists to collect r1 and r2 files
    r1s = []
    r2s = []

    # list the accession files for this sample
    fq_list = glob.glob('results/main/{cline}/reads/split_fastqs/*fastq'.format(cline=wildcards.cline))

    # parse through the accession lists and get the r1 and r2 paths
    for fq in fq_list:
        read_num = os.path.basename(fq).split('_')[-1].split('.fastq')[0]
        if read_num == '1':
            r1s.append(fq)
        else:
            r2s.append(fq)
    # return a dict to unpack
    d = {'r1s': r1s, 'r2s': r2s}
    return(d)



# Renaming because _1 and _2 in file names can caused a problem
# which meant I reverted to using _R1 and _R2 in the HiC-Pro configuration file.
# local rule
rule rename_before_hicpro_with_parallel: 
    input:
        unpack(get_split_r1_r2_fastqs),
        split_complete = 'results/main/{cline}/reads/split_fastqs/split.complete'
    output:
        new_dir = directory('results/main/{cline}/reads/renamed_fastqs_with_parallel/'),
        rename_complete = touch('results/main/{cline}/reads/renamed_fastqs_with_parallel/renamed.complete')
    log:
        'results/main/{cline}/logs/rule_rename_before_hicpro_with_parallel_{cline}.log'
    shell:
        r"""
            mkdir -p {output.new_dir}

            # renaming R1's
            echo "# renaming R1's" > {log}
            for fn in {input.r1s};
            do
                echo "fn: $fn"

                # get the new fn
                echo "	# get the new fn" >> {log} 
                new_fn=$(basename $fn | sed "s/_1\.fastq/_R1.fastq/")

                # get the current srr
                echo "	# get the current srr" >> {log}
                srr=$(basename $fn | sed "s/^.*_\([A-Za-z0-9-]*\)_.*$/\1/")

                # make an srr based directory
                echo "	# make an srr based directory" >> {log}
                new_dir="{output.new_dir}/${{srr}}"
                mkdir -p $new_dir

                # symlink the renamed fn to an srr based directory
                echo "	# symlink the renamed fn to an srr based directory" >> {log}
                new_fn="{output.new_dir}/${{srr}}/${{new_fn}}"
                abs_orig=$(readlink -f $fn)
                ln -s $abs_orig $new_fn 2>> {log}
            done
            
            # renaming R2's
            echo >> {log}
            echo "# renaming R2's" >> {log}
            for fn in {input.r2s};
            do
                echo "fn: $fn" >> {log}

                # get the new fn
                echo "	# get the new fn" >> {log}
                new_fn=$(basename $fn | sed "s/_2\.fastq/_R2.fastq/")

                # get the current srr
                echo "	# get the current srr" >> {log}
                srr=$(basename $fn | sed "s/^.*_\([A-Za-z0-9-]*\)_.*$/\1/")

                # symlink the renamed fn to an srr based directory
                echo "	# symlink the renamed fn to an srr based directory" >> {log}
                new_fn="{output.new_dir}/${{srr}}/${{new_fn}}"
                abs_orig=$(readlink -f $fn)

                echo $abs_orig >> {log}
                echo $new_fn >> {log}
                ln -s $abs_orig $new_fn 2>> {log}
            done
        """

# Align the HiC data with merging capability
# conda environments not working(?), left it for now
rule hicpro_with_parallel: # localrule
    input:
        fastq_dir = rules.rename_before_hicpro_with_parallel.output.new_dir,
        gs = rules.download_hg38_files.output.genome_sizes,
        digestion = re_digestion_file,
        bowtie2_idxs = rules.bowtie2_index_ref_genome.output,
        config = ancient(re_config_file),
    output:
        qsubs_written = touch('results/main/{cline}/hicpro_with_parallel/qsubs.written'),
    params:
        datadir = 'results/main/{cline}/hicpro/renamed_fastqs_with_parallel/', # part of rule rename_before_hicpr
        outdir = 'results/main/{cline}/hicpro_with_parallel/',
    log:
        'results/main/{cline}/logs/rule_hicpro_with_parallel_{cline}.log'
    conda:
        'envs/HiC-Pro-3.0.0.yml'
    shell:
        """
            # remove the snakemake made outdir since HiCPro wants to make it
            rm -r {params.outdir}

            # getting absoluate paths for data and outdirs, required
            # by HiCPro
            abs_datadir=$(readlink -f {input.fastq_dir})
            abs_outdir=$(readlink -f {params.outdir})

            # running with setting -s so that it runs only the alignment
            # pipeline (default settings), qsub jobs are automatically submitted
            echo "# running with setting -s so that it runs only the alignment" >> {log} 2>&1
            /home/jreyna/software/HiC-Pro-Full/HiC-Pro-3.0.0/bin/HiC-Pro \
                    -p \
                    -i $abs_datadir \
                    -o $abs_outdir \
                    -c {input.config} >> {log} 2>&1
        """


# After HiCPro genereates qsubs I submit them in a separate rule,
# for technical reasons I wasn't able to merge this rule
# and the preceeding rule, there seems to be some error signaling by
# HiC-Pro after it's done running which immediately stops the snakemake
# rule and nothing else can be run so this extra rule is essential. Also,
# I don't do any logging here because you have to cd into a different
# directory and that affects the pathing on the logs; it's such a simple
# rule that I don't think a log is necessary, you can make it work but 
# not work the effort.
rule hicpro_with_parallel_started: # localrule
    input:
        qsubs_written = rules.hicpro_with_parallel.output.qsubs_written
    output:
        qsubs_started = touch('results/main/{cline}/hicpro_with_parallel/qsubs.started')
    params:
        outdir = 'results/main/{cline}/hicpro_with_parallel/',
    conda:
        'envs/HiC-Pro-3.0.0.yml'
    shell:
        """
            # submit step 1
            echo "# submit step 1"
            cd {params.outdir} # need to run from the outdir, required by HiCPro
            qids=$(qsub HiCPro_step1_.sh)

            # submit step 2 with a hold
            echo "# submit step 2 with a hold"
            qsub -W depend=afterokarray:$qids HiCPro_step2_.sh
        """

