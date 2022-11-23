min_qual_iGenVar = list(range(config["quality_ranges"]["iGenVar"]["from"],
                              config["quality_ranges"]["iGenVar"]["to"],
                              config["quality_ranges"]["iGenVar"]["step"]))

rule filter_vcf:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.vcf"
    output:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'QUAL>={wildcards.min_qual}.00' {input.vcf} > {output.vcf}"

rule bgzip:
    input:
        "{name}.vcf"
    output:
        "{name}.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

rule tabix:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

rule truvari:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_iGenVar_only/{input_combination}/variants.min_qual_{min_qual}.vcf.gz.tbi"
    output:
        summary = "results/caller_comparison_iGenVar_only/eval/{input_combination}/min_qual_{min_qual}/summary.txt"
    log:
        "logs/truvari/truvari_output.{input_combination}.{min_qual}.log"
    params:
        output_dir = "results/caller_comparison_iGenVar_only/eval/{input_combination}/min_qual_{min_qual}"
    run:
        if wildcards.input_combination == 'S4': # NA12878 WGS GRCh38
            truth_set_gz = config["truth_set_NA12878"]["gz"],
            truth_set_bed = config["truth_set_NA12878"]["bed"]
        else:
            truth_set_gz = config["truth_set_HG002_renamed_chr"]["gz"],
            truth_set_bed = config["truth_set_HG002_renamed_chr"]["bed"]
        shell("""
            rm -rf {params.output_dir} && truvari bench -b {truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
                --passonly --includebed {truth_set_bed} &>> {log}
        """)

rule reformat_truvari_results:
    input:
        "results/caller_comparison_iGenVar_only/eval/{input_combination}/min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_iGenVar_only/eval/{input_combination}/min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.input_combination}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule cat_truvari_results_all:
    input:
        S1   = expand("results/caller_comparison_iGenVar_only/eval/S1/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2   = expand("results/caller_comparison_iGenVar_only/eval/S2/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S3   = expand("results/caller_comparison_iGenVar_only/eval/S3/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S4   = expand("results/caller_comparison_iGenVar_only/eval/S4/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S1L1 = expand("results/caller_comparison_iGenVar_only/eval/S1L1/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2L1 = expand("results/caller_comparison_iGenVar_only/eval/S2L1/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S3L1 = expand("results/caller_comparison_iGenVar_only/eval/S3L1/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S1L2 = expand("results/caller_comparison_iGenVar_only/eval/S1L2/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2L2 = expand("results/caller_comparison_iGenVar_only/eval/S2L2/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S3L2 = expand("results/caller_comparison_iGenVar_only/eval/S3L2/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S1L3 = expand("results/caller_comparison_iGenVar_only/eval/S1L3/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2L3 = expand("results/caller_comparison_iGenVar_only/eval/S2L3/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S3L3 = expand("results/caller_comparison_iGenVar_only/eval/S3L3/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L1   = expand("results/caller_comparison_iGenVar_only/eval/L1/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L2   = expand("results/caller_comparison_iGenVar_only/eval/L2/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L3   = expand("results/caller_comparison_iGenVar_only/eval/L3/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar)
    output:
        S1   = temp("results/caller_comparison_iGenVar_only/eval/S1.all_results.txt"),
        S2   = temp("results/caller_comparison_iGenVar_only/eval/S2.all_results.txt"),
        S3   = temp("results/caller_comparison_iGenVar_only/eval/S3.all_results.txt"),
        S4   = temp("results/caller_comparison_iGenVar_only/eval/S4.all_results.txt"),
        S1L1 = temp("results/caller_comparison_iGenVar_only/eval/S1L1.all_results.txt"),
        S2L1 = temp("results/caller_comparison_iGenVar_only/eval/S2L1.all_results.txt"),
        S3L1 = temp("results/caller_comparison_iGenVar_only/eval/S3L1.all_results.txt"),
        S1L2 = temp("results/caller_comparison_iGenVar_only/eval/S1L2.all_results.txt"),
        S2L2 = temp("results/caller_comparison_iGenVar_only/eval/S2L2.all_results.txt"),
        S3L2 = temp("results/caller_comparison_iGenVar_only/eval/S3L2.all_results.txt"),
        S1L3 = temp("results/caller_comparison_iGenVar_only/eval/S1L3.all_results.txt"),
        S2L3 = temp("results/caller_comparison_iGenVar_only/eval/S2L3.all_results.txt"),
        S3L3 = temp("results/caller_comparison_iGenVar_only/eval/S3L3.all_results.txt"),
        L1   = temp("results/caller_comparison_iGenVar_only/eval/L1.all_results.txt"),
        L2   = temp("results/caller_comparison_iGenVar_only/eval/L2.all_results.txt"),
        L3   = temp("results/caller_comparison_iGenVar_only/eval/L3.all_results.txt"),
        all = "results/caller_comparison_iGenVar_only/eval/all_results.txt"
    threads: 1
    run:
        shell("cat {input.S1} > {output.S1}")
        shell("cat {input.S2} > {output.S2}")
        shell("cat {input.S3} > {output.S3}")
        shell("cat {input.S4} > {output.S4}")
        shell("cat {input.S1L1} > {output.S1L1}")
        shell("cat {input.S2L1} > {output.S2L1}")
        shell("cat {input.S3L1} > {output.S3L1}")
        shell("cat {input.S1L2} > {output.S1L2}")
        shell("cat {input.S2L2} > {output.S2L2}")
        shell("cat {input.S3L2} > {output.S3L2}")
        shell("cat {input.S1L3} > {output.S1L3}")
        shell("cat {input.S2L3} > {output.S2L3}")
        shell("cat {input.S3L3} > {output.S3L3}")
        shell("cat {input.L1} > {output.L1}")
        shell("cat {input.L2} > {output.L2}")
        shell("cat {input.L3} > {output.L3}")
        shell("""
            cat {output.S1} {output.S2} {output.S3} {output.S4} \
                {output.S1L1} {output.S2L1} {output.S3L1} \
                {output.S1L2} {output.S2L2} {output.S3L2} \
                {output.S1L3} {output.S2L3} {output.S3L3} \
                {output.L1} {output.L2} {output.L3} > {output.all}
        """)
