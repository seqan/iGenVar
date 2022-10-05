sample = config["parameters"]["sample_HG002"],
sample_NA12878 = config["parameters"]["sample_NA12878"],
min_var_length = config["parameters"]["min_var_length"],
max_var_length = config["parameters"]["max_var_length"]

rule copy_iGenVar_results:
    output:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement,no_enhancement|full_enhancement|L2}/variants.vcf",
    run:
        if wildcards.long_read_enhancement == 'no_enhancement':
            if wildcards.dataset == 'Illumina_Paired_End':
                shell("cp results/caller_comparison_iGenVar_only/S1/variants.vcf {output.vcf}")
            else:
                shell("cp results/caller_comparison_iGenVar_only/S2/variants.vcf {output.vcf}")
        elif wildcards.long_read_enhancement == 'L2':
            shell("cp results/caller_comparison_iGenVar_only/L2/variants.vcf {output.vcf}") # PacBio CCS 30x coverage
        else: # wildcards.long_read_enhancement == 'full_enhancement':
            if wildcards.dataset == 'Illumina_Paired_End':
                shell("cp results/caller_comparison_iGenVar_only/S1L2/variants.vcf {output.vcf}")
            else:
                shell("cp results/caller_comparison_iGenVar_only/S2L2/variants.vcf {output.vcf}")
        # else: # L2x1 & L2x2 & L2x3 & ...

rule run_iGenVar:
    output:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement,SL2x1|SL2x2|SL2x3|SL2x5|SL2x10}/variants.vcf",
    log:
        "logs/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement,SL2x1|SL2x2|SL2x3|SL2x5|SL2x10}_output.log"
    threads: 2
    run:
        if wildcards.dataset == 'Illumina_Paired_End':
            short_bam = config["short_read_bam"]["s1"]
            # L2: The library was sequenced to approximately 30-fold coverage.
            if wildcards.long_read_enhancement == 'SL2x1':   # PacBio CCS 1x coverage
                long_bam = config["sampled_long_read_bam"]["l2x1"]
            elif wildcards.long_read_enhancement == 'SL2x2': # PacBio CCS 2x coverage
                long_bam = config["sampled_long_read_bam"]["l2x2"]
            elif wildcards.long_read_enhancement == 'SL2x3': # PacBio CCS 3x coverage
                long_bam = config["sampled_long_read_bam"]["l2x3"]
            elif wildcards.long_read_enhancement == 'SL2x5': # PacBio CCS 5x coverage
                long_bam = config["sampled_long_read_bam"]["l2x5"]
            else:        # long_read_enhancement == 'SL2x10': # PacBio CCS 10x coverage
                long_bam = config["sampled_long_read_bam"]["l2x10"]
            # else: # no_enhancement & L2 & ...
        else: # dataset == 'Illumina_Mate_Pair'
            short_bam = config["short_read_bam"]["s2"]
            if wildcards.long_read_enhancement == 'L2x1':   # PacBio CCS 1x coverage
                long_bam = config["sampled_long_read_bam"]["l2x1"]
            elif wildcards.long_read_enhancement == 'L2x2': # PacBio CCS 2x coverage
                long_bam = config["sampled_long_read_bam"]["l2x2"]
            elif wildcards.long_read_enhancement == 'L2x3': # PacBio CCS 3x coverage
                long_bam = config["sampled_long_read_bam"]["l2x3"]
            elif wildcards.long_read_enhancement == 'L2x5': # PacBio CCS 5x coverage
                long_bam = config["sampled_long_read_bam"]["l2x5"]
            else:        # long_read_enhancement == 'L2x10': # PacBio CCS 10x coverage
                long_bam = config["sampled_long_read_bam"]["l2x10"]
            # else: # no_enhancement & L2 & ...
        shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} --input_long_reads {long_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} --verbose \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
        """)
        # Defaults:
        # --method cigar_string --method split_read --method read_pairs --method read_depth
        # --clustering_methods hierarchical_clustering --refinement_methods no_refinement
        # --max_tol_inserted_length 50 --max_overlap 10 --hierarchical_clustering_cutoff 0.5

rule run_iGenVar_on_samples_only:
    output:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement,L2x1|L2x2|L2x3|L2x5|L2x10}/variants.vcf",
    log:
        "logs/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement,L2x1|L2x2|L2x3|L2x5|L2x10}_output.log"
    threads: 2
    run:
        # L2: The library was sequenced to approximately 30-fold coverage.
        if wildcards.long_read_enhancement == 'L2x1':   # PacBio CCS 1x coverage
            long_bam = config["sampled_long_read_bam"]["l2x1"]
        elif wildcards.long_read_enhancement == 'L2x2': # PacBio CCS 2x coverage
            long_bam = config["sampled_long_read_bam"]["l2x2"]
        elif wildcards.long_read_enhancement == 'L2x3': # PacBio CCS 3x coverage
            long_bam = config["sampled_long_read_bam"]["l2x3"]
        elif wildcards.long_read_enhancement == 'L2x5': # PacBio CCS 5x coverage
            long_bam = config["sampled_long_read_bam"]["l2x5"]
        else:        # long_read_enhancement == 'L2x10': # PacBio CCS 10x coverage
            long_bam = config["sampled_long_read_bam"]["l2x10"]
        shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_long_reads {long_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} --verbose \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
        """)
