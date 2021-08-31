rule run_igenvar:
    input:
        bam = config["long_bam"]
    output:
        vcf = "results/caller_comparison/iGenVar/variants.vcf"
    params:
        sample = config["parameters"]["sample"],
        min_qual = config["parameters"]["min_qual"],
        min_var_length = config["parameters"]["min_var_length"],
        max_var_length = config["parameters"]["max_var_length"]
    shell:
        """
        ./build/iGenVar/bin/iGenVar -t 1 -j {input.bam} -o {output.vcf} \
        --vcf_sample_name {params.sample} \
        --method cigar_string \
        --method split_read \
        --min_var_length {params.min_var_length} \
        --max_var_length {params.max_var_length} \
        --min_qual {params.min_qual}
        """
        # Defaults:
        # --clustering_methods hierarchical_clustering --refinement_methods no_refinement
        # --max_tol_inserted_length 5 --max_overlap 10 --hierarchical_clustering_cutoff 10

# SVIM
rule run_svim:
    input:
        bam = config["long_bam"],
        bai = config["long_bai"],
        genome = config["reference_fa_gz"]
    output:
        "results/caller_comparison/SVIM/variants.vcf"
    resources:
        mem_mb = 20000,
        time_min = 600,
        io_gb = 100
    params:
        working_dir = "results/caller_comparison/SVIM/",
        sample = config["parameters"]["sample"],
        min_var_length = config["parameters"]["min_var_length"],
        max_var_length = config["parameters"]["max_var_length"]
    threads: 1
    conda:
        "../../../envs/svim.yaml"
    shell:
        """
        svim alignment --sample {params.sample} \
        --partition_max_distance 1000 \
        --cluster_max_distance 0.5 \
        --min_sv_size {params.min_var_length} \
        --segment_gap_tolerance 20 \
        --segment_overlap_tolerance 20 \
        --interspersed_duplications_as_insertions \
        --tandem_duplications_as_insertions \
        --read_names \
        --max_sv_size {params.max_var_length} \
        --verbose \
        {params.working_dir} {input.bam} {input.genome}
        """
        # Defaults:
        # --position_distance_normalizer 900 --edit_distance_normalizer 1.0

# SNIFFLES (we have to loop over min_support, because sniffles does not write a quality score into the vcf)
rule run_sniffles:
    input:
        bam = config["long_bam_md"],
    output:
        expand("results/caller_comparison/Sniffles/raw_variants.{minsupport}.vcf",
               minsupport=list(range(config["minimums"]["sniffles_from"],
                                     config["minimums"]["sniffles_to"]+1,
                                     config["minimums"]["sniffles_step"])))
    resources:
        mem_mb = 400000,
        time_min = 1200,
        io_gb = 100
    params:
        min_support = config["parameters"]["min_qual"],
        min_length = config["parameters"]["min_var_length"],
        qual_from = config["minimums"]["sniffles_from"],
        qual_to = config["minimums"]["sniffles_to"]+1,
        qual_step = config["minimums"]["sniffles_step"]
    threads: 10
    conda:
        "../../../envs/sniffles.yaml"
    shell:
        """
        for i in $(seq {params.qual_from} {params.qual_step} {params.qual_to})
        do
            sniffles --mapped_reads {input.bam} --vcf results/caller_comparison/Sniffles/raw_variants.$i.vcf \
            --min_support $i --min_length {params.min_length} --threads {threads} --genotype
        done
        """

#see https://github.com/spiralgenetics/truvari/issues/43
rule fix_sniffles:
    input:
        "results/caller_comparison/Sniffles/raw_variants.{support}.vcf"
    output:
        "results/caller_comparison/Sniffles/variants.unsorted.min_qual_{support,[0-9]+}.vcf"
    shell:
        "sed 's/##INFO=<ID=SUPTYPE,Number=A/##INFO=<ID=SUPTYPE,Number=./' {input} > {output}"

# Split to SV classes
# Since iGenVar can only find INS and DEL so far, we filter these out for better comparability.
rule fix_sniffles_2_and_filter_insertions_and_deletions:
    input:
        "results/caller_comparison/Sniffles/variants.unsorted.min_qual_{support,[0-9]+}.vcf"
    output:
        "results/caller_comparison/Sniffles/variants.min_qual_{support,[0-9]+}.vcf"
    shell:
        "bcftools view -i 'SVTYPE=\"DEL\" | SVTYPE=\"INS\"' {input} | bcftools sort > {output}"

#PBSV
rule run_pbsv_dicsover:
    input:
        bam = config["long_bam"]
    output:
        svsig_gz = dynamic("results/caller_comparison/pbsv/signatures.{region}.svsig.gz")
    resources:
        mem_mb = 400000,
        time_min = 2000,
        io_gb = 100
    threads: 1
    conda:
        "../../../envs/pbsv.yaml"
    shell:
        # "pbsv discover {input.bam} {output.svsig_gz}"
        """
        for i in $(samtools view -H {input.bam} | grep '^@SQ' | cut -f2 | cut -d':' -f2); do
            pbsv discover --region $i {input.bam} results/caller_comparison/pbsv/signatures.$i.svsig.gz
        done
        """

rule run_pbsv_call:
    input:
        genome = config["reference_fa"],
        svsig_gz = dynamic("results/caller_comparison/pbsv/signatures.{region}.svsig.gz")
    output:
        vcf = expand("results/caller_comparison/pbsv/variants.min_qual_{minsupport}.vcf",
                     minsupport=list(range(config["minimums"]["pbsv_from"],
                                           config["minimums"]["pbsv_to"]+1,
                                           config["minimums"]["pbsv_step"])))
    resources:
        mem_mb = 400000,
        time_min = 2000,
        io_gb = 100
    params:
        min_sv_length = config["parameters"]["min_var_length"],
        qual_from = config["minimums"]["pbsv_from"],
        qual_to = config["minimums"]["pbsv_to"]+1,
        qual_step = config["minimums"]["pbsv_step"]
    threads: 1
    conda:
        "../../../envs/pbsv.yaml"
    shell:
        # pbsv call --types DEL,INS,DUP --min-sv-length {params.min_sv_length} --max-ins-length 100K \
        """
        for i in $(seq {params.qual_from} {params.qual_step} {params.qual_to})
        do
            pbsv call --types DEL,INS --min-sv-length {params.min_sv_length} --max-ins-length 100K \
            --call-min-reads-all-samples $i --call-min-reads-one-sample $i \
            --call-min-reads-per-strand-all-samples 0 --call-min-bnd-reads-all-samples 0 --call-min-read-perc-one-sample 0 \
            --num-threads {threads} {input.genome} {input.svsig_gz} results/caller_comparison/pbsv/variants.min_qual_$i.vcf
        done
        """
