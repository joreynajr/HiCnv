# Renaming because _1 and _2 in file names can caused a problem
# which meant I reverted to using _R1 and _R2 in the HiC-Pro configuration file.
rule rename_before_hicpro_with_parallel:
    input:
        unpack(get_r1_r2_fastqs)
    output:
        new_dir = directory('results/main/{cline}/reads/renamed_fastqs_with_parallel/'),
        rename_complete = touch('results/main/{cline}/reads/renamed_fastqs_with_parallel/renamed.complete')
    log:
        'results/main/{cline}/logs/rule_rename_before_hicpro_with_parallel_{cline}.log'
    shell:
        r"""
            mkdir -p {output.new_dir}

            # renaming R1's
            for fn in {input.r1s};
            do
                # get the new fn
                new_fn=$(basename $fn | sed "s/_1\.fastq\.gz/_R1.fastq.gz/")

                # get the current srr
                srr=$(basename $fn | sed "s/_.*//")

                # make an srr based directory
                new_dir="{output.new_dir}/${{srr}}"
                mkdir $new_dir

                # symlink the renamed fn to an srr based directory
                new_fn="{output.new_dir}/${{srr}}/${{new_fn}}"
                abs_orig=$(readlink -f $fn)
                ln -s $abs_orig $new_fn
            done

            # renaming R2's
            for fn in {input.r2s};
            do
                # get the new fn
                new_fn=$(basename $fn | sed "s/_2\.fastq\.gz/_R2.fastq.gz/")

                # get the current srr
                srr=$(basename $fn | sed "s/_.*//")

                # symlink the renamed fn to an srr based directory
                new_fn="{output.new_dir}/${{srr}}/${{new_fn}}"
                abs_orig=$(readlink -f $fn)
                ln -s $abs_orig $new_fn
            done
        """


# Align the HiC data with merging capability
# conda environments not working(?), left it for now
rule hicpro_align_only_with_parallel_all: # localrule
    input:
        fastq_dir = rules.rename_before_hicpro_with_parallel.output.new_dir,
        gs = rules.download_hg38_files.output.genome_sizes,
        digestion = re_digestion_file,
        bowtie2_idxs = rules.bowtie2_index_ref_genome.output,
        config = ancient(re_config_file),
    output:
        bowtie_running = touch('results/main/{cline}/hicpro_with_parallel/bowtie_results/all.running')
    params:
        datadir = 'results/main/{cline}/hicpro/renamed_fastqs_with_parallel/', # part of rule rename_before_hicpr
        outdir = 'results/main/{cline}/hicpro_with_parallel/',
    log:
        'results/main/{cline}/logs/rule_hicpro_align_with_parallel_only_{cline}.log'
    benchmark:
        'results/main/{cline}/benchmarks/rule_hicpro_align_only_with_parallel_{cline}.bmk'
    conda:
        'workflow/envs/HiC-Pro-3.0.0.yml'
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
                    -c {input.config} >> {log} 2>&1 || cd {params.outdir}

            # submit step 1 # automatic submissions are not working, will improve in the future
            #echo "# submit step 1" >> {log} 2>&1
            #qids=$(qsub -w {params.outdir} {params.outdir}/HiCPro_step1_.sh) >> {log} 2>&1

            # submit step 2 with a hold; same problem as qsub with step 1
            #echo "# submit step 2 with a hold" >> {log} 2>&1
            #qsub -w {params.outdir} -W depend=afterokarray:$qids {params.outdir}/HiCPro_step2_.sh >> {log} 2>&1
        """


# Splitting the original R1 and R2 fastq's because they are really big
# Current this is set up to run serially but by adding an srr wildcard
# I can set it up to run in a more parallel way. The idea right now
# is that you would have a single large >70gb file per sample (bio-rep).
rule split_before_hicpro:
    input:
        unpack(get_r1_r2_fastqs)
    output:
        outdir = directory('results/main/{cline}/reads/split_fastqs/'),
        split_complete = touch('results/main/{cline}/reads/split_fastqs/split.complete')
    params:
        nreads = 50000000
    log:
        'results/main/{cline}/logs/rule_split_before_hicpro_{cline}.log'
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
