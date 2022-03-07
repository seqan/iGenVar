sample = config["parameters"]["sample"],
min_var_length = config["parameters"]["min_var_length"],
max_var_length = config["parameters"]["max_var_length"]

# iGenVar
# TODO (irallia 07.03.2022): change this to short read input
rule run_iGenVar_S:
    input:
        bam = config["long_bam"]
    output:
        vcf = "results/caller_comparison_short_read/iGenVar/variants.vcf"
    threads: 1
    shell:
        """
        ./build/iGenVar/bin/iGenVar --input_long_reads {input.bam} --output {output.vcf} --vcf_sample_name {sample} \
            --threads {threads} --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1
        """

# TODO (irallia 07.03.2022): add short read input
rule run_iGenVar_SL:
    input:
        bam = config["long_bam"]
    output:
        vcf = "results/caller_comparison_short_read/iGenVar/variants.vcf"
    threads: 1
    shell:
        """
        ./build/iGenVar/bin/iGenVar --input_long_reads {input.bam} --output {output.vcf} --vcf_sample_name {sample} \
            --threads {threads} --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1
        """
        # Defaults:
        # --method cigar_string --method split_read --method read_pairs --method read_depth
        # --clustering_methods hierarchical_clustering --refinement_methods no_refinement
        # --max_tol_inserted_length 50 --max_overlap 10 --hierarchical_clustering_cutoff 0.5
