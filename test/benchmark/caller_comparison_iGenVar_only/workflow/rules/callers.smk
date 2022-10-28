sample = config["parameters"]["sample"],
min_var_length = config["parameters"]["min_var_length"],
max_var_length = config["parameters"]["max_var_length"]

rule run_iGenVar:
    output:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/unsorted_variants.vcf"
    log:
        "logs/caller_comparison_iGenVar_only/{input_combination}_output.log"
    threads: 4
    run:
        if wildcards.input_combination == 'S1': # Illumina Paired End
            short_bam = config["short_read_bam"]["s1"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'S2': # Illumina Mate Pair
            short_bam = config["short_read_bam"]["s2"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'S1L1': # Illumina Paired End & MtSinai PacBio
            short_bam = config["short_read_bam"]["s1"],
            long_bam = config["long_read_bam"]["l1"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} --input_long_reads {long_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'S2L1': # Illumina Mate Pair & MtSinai PacBio
            short_bam = config["short_read_bam"]["s2"],
            long_bam = config["long_read_bam"]["l1"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} --input_long_reads {long_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'S1L2': # Illumina Paired End & PacBio CCS
            short_bam = config["short_read_bam"]["s1"],
            long_bam = config["long_read_bam"]["l2"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} --input_long_reads {long_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'S2L2': # Illumina Mate Pair & PacBio CCS
            short_bam = config["short_read_bam"]["s2"],
            long_bam = config["long_read_bam"]["l2"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} --input_long_reads {long_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'S1L3': # Illumina Paired End & 10X Genomics
            short_bam = config["short_read_bam"]["s1"],
            long_bam = config["long_read_bam"]["l3"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} --input_long_reads {long_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'S2L3': # Illumina Mate Pair & 10X Genomics
            short_bam = config["short_read_bam"]["s2"],
            long_bam = config["long_read_bam"]["l3"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} --input_long_reads {long_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'L1': # MtSinai PacBio
            long_bam = config["long_read_bam"]["l1"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_long_reads {long_bam} \
                 --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'L2': # PacBio CCS
            long_bam = config["long_read_bam"]["l2"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_long_reads {long_bam} \
                 --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'L3': # 10X Genomics
            long_bam = config["long_read_bam"]["l3"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_long_reads {long_bam} \
                 --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'hg38_Sim_default': # Illumina
            short_bam = config["simulated_short_read_bam"]["sim1"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'hg38_Sim_InDel': # Illumina
            short_bam = config["simulated_short_read_bam"]["sim2"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        elif wildcards.input_combination == 'hg38_Sim_noSNP': # Illumina
            short_bam = config["simulated_short_read_bam"]["sim3"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} \
                --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        else: # wildcards.input_combination == 'hg38_Sim_SNPandSV': # Illumina
            short_bam = config["simulated_short_read_bam"]["sim4"]
            shell("""
            /usr/bin/time -v ./build/iGenVar/bin/iGenVar --input_short_reads {short_bam} \
                 --output {output.vcf} --vcf_sample_name {sample} --threads {threads} \
                --min_var_length {min_var_length} --max_var_length {max_var_length} --min_qual 1 &>> {log}
            """)
        # Defaults:
        # --method cigar_string --method split_read --method read_pairs --method read_depth
        # --clustering_methods hierarchical_clustering --refinement_methods no_refinement
        # --max_tol_inserted_length 50 --max_overlap 10 --hierarchical_clustering_cutoff 0.5

# TODO (irallia 24.06.2022): The vcf is not sorted anymore.
rule picard:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/unsorted_variants.vcf"
    output:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.vcf"
    log:
        "logs/caller_comparison_iGenVar_only/picard_output.{input_combination}.log"
    conda:
        "../../../envs/simulation.yaml"
    shell:
        "picard SortVcf -I {input.vcf} -O {output.vcf} -Xms1g -Xmx100g --TMP_DIR tmp/picard/ &>> {log}"
