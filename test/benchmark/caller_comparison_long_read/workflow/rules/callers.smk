sample = config["parameters"]["sample"],
min_var_length = config["parameters"]["min_var_length"],
max_var_length = config["parameters"]["max_var_length"]

rule run_igenvar:
    input:
        bam = config["long_read_bam"]["l2"]
    output:
        vcf = "results/caller_comparison_long_read/iGenVar/variants.vcf"
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

# SVIM
rule run_svim:
    input:
        bam = config["long_read_bam"]["l2"],
        genome = config["reference_fa"]["PacBio_CCS"] # [E::fai_build3_core] Cannot index files compressed with gzip, please use bgzip
    output:
        "results/caller_comparison_long_read/SVIM/variants.vcf"
    resources:
        mem_mb = 20000,
        time_min = 600,
        io_gb = 100
    params:
        working_dir = "results/caller_comparison_long_read/SVIM/"
    threads: 1
    # conda:
    #     "../../../envs/svim.yaml"
    shell:
        """
        svim alignment --sample {sample} --partition_max_distance 1000 --cluster_max_distance 0.5 \
            --min_sv_size {min_var_length} --segment_gap_tolerance 20 --segment_overlap_tolerance 20 --read_names \
            --max_sv_size {max_var_length} {params.working_dir} {input.bam} {input.genome} &>> logs/svim_output.log
        """
        # Defaults:
        # --position_distance_normalizer 900 --edit_distance_normalizer 1.0

# SNIFFLES (we have to loop over min_support, because sniffles does not write a quality score into the vcf)
rule run_sniffles:
    input:
        bam = config["long_md_bam"]
    output:
        "results/caller_comparison_long_read/Sniffles/raw_variants_{min_qual}.vcf"
    resources:
        mem_mb = 400000,
        time_min = 1200,
        io_gb = 100
    threads: 1
    conda:
        "../../../envs/sniffles.yaml"
    shell:
        """
        sniffles --mapped_reads {input.bam} \
            --vcf results/caller_comparison_long_read/Sniffles/raw_variants_{wildcards.min_qual}.vcf \
            --min_support {wildcards.min_qual} --min_length {min_var_length} --threads {threads} --genotype \
            >> logs/sniffles_output.log
        """

# see https://github.com/spiralgenetics/truvari/issues/43
# see https://github.com/fritzsedlazeck/Sniffles/issues/209 (Fixed, but there just hasn't been a release since the fix.)
rule fix_sniffles:
    input:
        "results/caller_comparison_long_read/Sniffles/raw_variants_{min_qual}.vcf"
    output:
        "results/caller_comparison_long_read/Sniffles/variants.unsorted.min_qual_{min_qual}.vcf"
    run:
        shell("sed 's/##INFO=<ID=SUPTYPE,Number=A/##INFO=<ID=SUPTYPE,Number=./' {input} > {output}")
        shell("sed -i '4i##FILTER=<ID=STRANDBIAS,Description=\"Strand is biased.\">' {output}")

# Split to SV classes
rule fix_sniffles_2:
    input:
        "results/caller_comparison_long_read/Sniffles/variants.unsorted.min_qual_{min_qual}.vcf"
    output:
        "results/caller_comparison_long_read/Sniffles/variants.min_qual_{min_qual}.vcf"
    shell:
        "bcftools sort {input} > {output}"

#PBSV
rule run_pbsv_dicsover:
    input:
        bam = config["long_read_bam"]["l2"]
    output:
        svsig_gz = "results/caller_comparison_long_read/pbsv/signatures.svsig.gz"
        # svsig_gz = dynamic("results/caller_comparison_long_read/pbsv/signatures.{region}.svsig.gz")
    resources:
        mem_mb = 400000,
        time_min = 2000,
        io_gb = 100
    threads: 1
    conda:
        "../../../envs/pbsv.yaml"
    shell:
        "pbsv discover {input.bam} {output.svsig_gz}"
        # """
        # for i in $(samtools view -H {input.bam} | grep '^@SQ' | cut -f2 | cut -d':' -f2); do
        #     pbsv discover --region $i {input.bam} results/caller_comparison_long_read/pbsv/signatures.$i.svsig.gz
        # done
        # """

rule run_pbsv_call:
    input:
        genome = config["reference_fa"]["PacBio_CCS"],
        svsig_gz = "results/caller_comparison_long_read/pbsv/signatures.svsig.gz"
        # svsig_gz = dynamic("results/caller_comparison_long_read/pbsv/signatures.{region}.svsig.gz")
    output:
        "results/caller_comparison_long_read/pbsv/variants.min_qual_{min_qual}.vcf"
    resources:
        mem_mb = 400000,
        time_min = 2000,
        io_gb = 100
    threads: 1
    conda:
        "../../../envs/pbsv.yaml"
    shell:
        """
        pbsv call --types DEL,INS,DUP --min-sv-length {min_var_length} --max-ins-length 100K \
            --call-min-reads-all-samples {wildcards.min_qual} --call-min-reads-one-sample {wildcards.min_qual} \
            --call-min-reads-per-strand-all-samples 0 --call-min-bnd-reads-all-samples 0 --call-min-read-perc-one-sample 0 \
            --num-threads {threads} {input.genome} {input.svsig_gz} \
            results/caller_comparison_long_read/pbsv/variants.min_qual_{wildcards.min_qual}.vcf
        """

rule run_pbsv_call_without_DUP:
    input:
        genome = config["reference_fa"]["PacBio_CCS"],
        svsig_gz = "results/caller_comparison_long_read/pbsv/signatures.svsig.gz"
        # svsig_gz = dynamic("results/caller_comparison_long_read/pbsv/signatures.{region}.svsig.gz")
    output:
        "results/caller_comparison_long_read/pbsv_without_DUP/variants.min_qual_{min_qual}.vcf"
    resources:
        mem_mb = 400000,
        time_min = 2000,
        io_gb = 100
    threads: 1
    conda:
        "../../../envs/pbsv.yaml"
    shell:
        """
        pbsv call --types DEL,INS --min-sv-length {min_var_length} --max-ins-length 100K \
            --call-min-reads-all-samples {wildcards.min_qual} --call-min-reads-one-sample {wildcards.min_qual} \
            --call-min-reads-per-strand-all-samples 0 --call-min-bnd-reads-all-samples 0 --call-min-read-perc-one-sample 0 \
            --num-threads {threads} {input.genome} {input.svsig_gz} \
            results/caller_comparison_long_read/pbsv_without_DUP/variants.min_qual_{wildcards.min_qual}.vcf
        """
