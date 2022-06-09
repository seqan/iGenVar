min_qual_VaquitaLR = list(range(config["quality_ranges"]["Vaquita-LR"]["from"],
                                config["quality_ranges"]["Vaquita-LR"]["to"],
                                config["quality_ranges"]["Vaquita-LR"]["step"]))

rule DUP_as_INS:
    input:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.vcf"
    output:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.DUP_as_INS.vcf"
    run:
        shell("sed -e 's/<DUP:TANDEM>/<INS>/g' {input.vcf} | sed -e 's/SVTYPE=DUP/SVTYPE=INS/g' > {output.vcf}")
        shell("sed -e 's/<DUP>/<INS>/g' {input.vcf} | sed -e 's/SVTYPE=DUP/SVTYPE=INS/g' > {output.vcf}")

rule DUP_as_INS_filter_vcf:
    input:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.DUP_as_INS.vcf"
    output:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.DUP_as_INS.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'INFO/SC>={wildcards.min_qual}' {input.vcf} > {output.vcf}"

rule DUP_as_INS_truvari:
    input:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.DUP_as_INS.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_vaquita_lr/{input_combination}/variants.DUP_as_INS.min_qual_{min_qual}.vcf.gz.tbi"
    output:
        summary = "results/caller_comparison_vaquita_lr/eval/{input_combination}/DUP_as_INS.min_qual_{min_qual}/summary.txt"
    log:
        "logs/truvari/truvari_output.{input_combination}.DUP_as_INS.{min_qual}.log"
    params:
        output_dir = "results/caller_comparison_vaquita_lr/eval/{input_combination}/DUP_as_INS.min_qual_{min_qual}"
    run:
        if (wildcards.input_combination == 'S2') | (wildcards.input_combination == 'L3') | (wildcards.input_combination == 'S2L3'):
            truth_set_gz = config["truth_set_renamed_chr"]["gz"],
            truth_set_bed = config["truth_set_renamed_chr"]["bed"]
        else: # S1, L1, L2, S1L1, S1L2
            truth_set_gz = config["truth_set"]["gz"],
            truth_set_bed = config["truth_set"]["bed"]
        shell("""
        rm -rf {params.output_dir} && truvari bench -b {truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
            --passonly --includebed {truth_set_bed} &>> {log}
        """)
        # => The combinations of S1L3, S2L1, S2L2 are not possible.
        # -f data/reference/hs37d5.fa

rule DUP_as_INS_reformat_truvari_results:
    input:
        "results/caller_comparison_vaquita_lr/eval/{input_combination}/DUP_as_INS.min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_vaquita_lr/eval/{input_combination}/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.input_combination}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule DUP_as_INS_cat_truvari_results_all:
    input:
        S1   = expand("results/caller_comparison_vaquita_lr/eval/S1/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        # S2   = expand("results/caller_comparison_vaquita_lr/eval/S2/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
        #               min_qual = min_qual_VaquitaLR),
        S1L1 = expand("results/caller_comparison_vaquita_lr/eval/S1L1/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        S1L2 = expand("results/caller_comparison_vaquita_lr/eval/S1L2/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        # S2L3 = expand("results/caller_comparison_vaquita_lr/eval/S2L3/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
        #               min_qual = min_qual_VaquitaLR),
        L1   = expand("results/caller_comparison_vaquita_lr/eval/L1/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        L2   = expand("results/caller_comparison_vaquita_lr/eval/L2/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        L3   = expand("results/caller_comparison_vaquita_lr/eval/L3/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR)
    output:
        S1   = temp("results/caller_comparison_vaquita_lr/eval/S1.DUP_as_INS.all_results.txt"),
        # S2   = temp("results/caller_comparison_vaquita_lr/eval/S2.DUP_as_INS.all_results.txt"),
        S1L1 = temp("results/caller_comparison_vaquita_lr/eval/S1L1.DUP_as_INS.all_results.txt"),
        S1L2 = temp("results/caller_comparison_vaquita_lr/eval/S1L2.DUP_as_INS.all_results.txt"),
        # S2L3 = temp("results/caller_comparison_vaquita_lr/eval/S2L3.DUP_as_INS.all_results.txt"),
        L1   = temp("results/caller_comparison_vaquita_lr/eval/L1.DUP_as_INS.all_results.txt"),
        L2   = temp("results/caller_comparison_vaquita_lr/eval/L2.DUP_as_INS.all_results.txt"),
        L3   = temp("results/caller_comparison_vaquita_lr/eval/L3.DUP_as_INS.all_results.txt"),
        all = "results/caller_comparison_vaquita_lr/eval/DUP_as_INS.all_results.txt"
    threads: 1
    run:
        shell("cat {input.S1} > {output.S1}")
        # shell("cat {input.S2} > {output.S2}")
        shell("cat {input.S1L1} > {output.S1L1}")
        shell("cat {input.S1L2} > {output.S1L2}")
        # shell("cat {input.S2L3} > {output.S2L3}")
        shell("cat {input.L1} > {output.L1}")
        shell("cat {input.L2} > {output.L2}")
        shell("cat {input.L3} > {output.L3}")
        shell("""
            cat {output.S1} \
                {output.S1L1} {output.S1L2} \
                {output.L1} {output.L2} {output.L3} > {output.all}
        """)
        # shell("""
        #     cat {output.S1} {output.S2} \
        #         {output.S1L1} {output.S1L2} {output.S2L3} \
        #         {output.L1} {output.L2} {output.L3} > {output.all}
        # """)
