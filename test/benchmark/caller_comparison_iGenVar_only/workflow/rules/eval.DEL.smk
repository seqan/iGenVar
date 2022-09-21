min_qual_iGenVar = list(range(config["quality_ranges"]["iGenVar"]["from"],
                              config["quality_ranges"]["iGenVar"]["to"],
                              config["quality_ranges"]["iGenVar"]["step"]))

rule filter_for_DEL:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.vcf"
    output:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.DEL.vcf"
    shell:
        "sed '/SVTYPE=DUP/d' {input.vcf} | sed '/SVTYPE=INS/d' | sed '/SVTYPE=INV/d' > {output.vcf}"

rule filter_DEL_vcf:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.DEL.vcf"
    output:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.DEL.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'QUAL>={wildcards.min_qual}.00' {input.vcf} > {output.vcf}"

rule truvari_DEL:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.DEL.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_iGenVar_only/{input_combination}/variants.DEL.min_qual_{min_qual}.vcf.gz.tbi"
    output:
        summary = "results/caller_comparison_iGenVar_only/eval/{input_combination}/DEL.min_qual_{min_qual}/summary.txt"
    log:
        "logs/truvari/truvari_output.{input_combination}.{min_qual}.log"
    params:
        output_dir = "results/caller_comparison_iGenVar_only/eval/{input_combination}/DEL.min_qual_{min_qual}",
        truth_set_gz = config["truth_set_DEL_renamed_chr"]["gz"],
        truth_set_bed = config["truth_set_renamed_chr"]["bed"]
    shell:
        """
        rm -rf {params.output_dir} && truvari bench -b {params.truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
            --passonly --includebed {params.truth_set_bed} &>> {log}
        """
        # -f data/reference/hs37d5.fa

rule reformat_truvari_results_DEL:
    input:
        "results/caller_comparison_iGenVar_only/eval/{input_combination}/DEL.min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_iGenVar_only/eval/{input_combination}/DEL.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.input_combination}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule cat_truvari_results_DEL:
    input:
        S1   = expand("results/caller_comparison_iGenVar_only/eval/S1/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2   = expand("results/caller_comparison_iGenVar_only/eval/S2/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S1L1 = expand("results/caller_comparison_iGenVar_only/eval/S1L1/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2L1 = expand("results/caller_comparison_iGenVar_only/eval/S2L1/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S1L2 = expand("results/caller_comparison_iGenVar_only/eval/S1L2/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2L2 = expand("results/caller_comparison_iGenVar_only/eval/S2L2/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S1L3 = expand("results/caller_comparison_iGenVar_only/eval/S1L3/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2L3 = expand("results/caller_comparison_iGenVar_only/eval/S2L3/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L1   = expand("results/caller_comparison_iGenVar_only/eval/L1/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L2   = expand("results/caller_comparison_iGenVar_only/eval/L2/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L3   = expand("results/caller_comparison_iGenVar_only/eval/L3/DEL.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar)
    output:
        S1   = temp("results/caller_comparison_iGenVar_only/eval/S1.DEL_results.txt"),
        S2   = temp("results/caller_comparison_iGenVar_only/eval/S2.DEL_results.txt"),
        S1L1 = temp("results/caller_comparison_iGenVar_only/eval/S1L1.DEL_results.txt"),
        S2L1 = temp("results/caller_comparison_iGenVar_only/eval/S2L1.DEL_results.txt"),
        S1L2 = temp("results/caller_comparison_iGenVar_only/eval/S1L2.DEL_results.txt"),
        S2L2 = temp("results/caller_comparison_iGenVar_only/eval/S2L2.DEL_results.txt"),
        S1L3 = temp("results/caller_comparison_iGenVar_only/eval/S1L3.DEL_results.txt"),
        S2L3 = temp("results/caller_comparison_iGenVar_only/eval/S2L3.DEL_results.txt"),
        L1   = temp("results/caller_comparison_iGenVar_only/eval/L1.DEL_results.txt"),
        L2   = temp("results/caller_comparison_iGenVar_only/eval/L2.DEL_results.txt"),
        L3   = temp("results/caller_comparison_iGenVar_only/eval/L3.DEL_results.txt"),
        all = "results/caller_comparison_iGenVar_only/eval/DEL_results.txt"
    threads: 1
    run:
        shell("cat {input.S1} > {output.S1}")
        shell("cat {input.S2} > {output.S2}")
        shell("cat {input.S1L1} > {output.S1L1}")
        shell("cat {input.S2L1} > {output.S2L1}")
        shell("cat {input.S1L2} > {output.S1L2}")
        shell("cat {input.S2L2} > {output.S2L2}")
        shell("cat {input.S1L3} > {output.S1L3}")
        shell("cat {input.S2L3} > {output.S2L3}")
        shell("cat {input.L1} > {output.L1}")
        shell("cat {input.L2} > {output.L2}")
        shell("cat {input.L3} > {output.L3}")
        shell("""
            cat {output.S1} {output.S2} {output.S1L1} {output.S2L1} \
                {output.S1L2} {output.S2L2} {output.S1L3} {output.S2L3} \
                {output.L1} {output.L2} {output.L3} > {output.all}
        """)
