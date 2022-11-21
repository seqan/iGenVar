min_qual_iGenVar = list(range(config["quality_ranges"]["iGenVar"]["from"],
                              config["quality_ranges"]["iGenVar"]["to"],
                              config["quality_ranges"]["iGenVar"]["step"]))

rule filter_for_INS:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.vcf"
    output:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.INS.vcf"
    shell:
        "sed '/SVTYPE=DEL/d' {input.vcf} | sed '/SVTYPE=DUP/d' | sed '/SVTYPE=INV/d' > {output.vcf}"

rule filter_INS_vcf:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.INS.vcf"
    output:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.INS.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'QUAL>={wildcards.min_qual}.00' {input.vcf} > {output.vcf}"

rule truvari_INS:
    input:
        vcf = "results/caller_comparison_iGenVar_only/{input_combination}/variants.INS.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_iGenVar_only/{input_combination}/variants.INS.min_qual_{min_qual}.vcf.gz.tbi"
    output:
        summary = "results/caller_comparison_iGenVar_only/eval/{input_combination}/INS.min_qual_{min_qual}/summary.txt"
    log:
        "logs/truvari/truvari_output.{input_combination}.{min_qual}.log"
    params:
        output_dir = "results/caller_comparison_iGenVar_only/eval/{input_combination}/INS.min_qual_{min_qual}",
        truth_set_gz = config["truth_set_INS_renamed_chr"]["gz"],
        truth_set_bed = config["truth_set_renamed_chr"]["bed"]
    shell:
        """
        rm -rf {params.output_dir} && truvari bench -b {params.truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
            --passonly --includebed {params.truth_set_bed} &>> {log}
        """
        # -f data/reference/hs37d5.fa

rule reformat_truvari_results_INS:
    input:
        "results/caller_comparison_iGenVar_only/eval/{input_combination}/INS.min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_iGenVar_only/eval/{input_combination}/INS.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.input_combination}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule cat_truvari_results_INS:
    input:
        S1   = expand("results/caller_comparison_iGenVar_only/eval/S1/INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        S2   = expand("results/caller_comparison_iGenVar_only/eval/S2/INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L1   = expand("results/caller_comparison_iGenVar_only/eval/L1/INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L2   = expand("results/caller_comparison_iGenVar_only/eval/L2/INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar),
        L3   = expand("results/caller_comparison_iGenVar_only/eval/L3/INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_iGenVar)
    output:
        S1   = temp("results/caller_comparison_iGenVar_only/eval/S1.INS_results.txt"),
        S2   = temp("results/caller_comparison_iGenVar_only/eval/S2.INS_results.txt"),
        L1   = temp("results/caller_comparison_iGenVar_only/eval/L1.INS_results.txt"),
        L2   = temp("results/caller_comparison_iGenVar_only/eval/L2.INS_results.txt"),
        L3   = temp("results/caller_comparison_iGenVar_only/eval/L3.INS_results.txt"),
        all = "results/caller_comparison_iGenVar_only/eval/INS_results.txt"
    threads: 1
    run:
        shell("cat {input.S1} > {output.S1}")
        shell("cat {input.S2} > {output.S2}")
        shell("cat {input.L1} > {output.L1}")
        shell("cat {input.L2} > {output.L2}")
        shell("cat {input.L3} > {output.L3}")
        shell("""
            cat {output.S1} {output.S2} \
                {output.L1} {output.L2} {output.L3} > {output.all}
        """)
