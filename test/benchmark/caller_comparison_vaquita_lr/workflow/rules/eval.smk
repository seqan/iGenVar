min_qual_VaquitaLR = list(range(config["quality_ranges"]["Vaquita-LR"]["from"],
                                config["quality_ranges"]["Vaquita-LR"]["to"],
                                config["quality_ranges"]["Vaquita-LR"]["step"]))

rule filter_vcf:
    input:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.vcf"
    output:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'INFO/SC>={wildcards.min_qual}' {input.vcf} > {output.vcf}"

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
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_vaquita_lr/{input_combination}/variants.min_qual_{min_qual}.vcf.gz.tbi"
    output:
        summary = "results/caller_comparison_vaquita_lr/eval/{input_combination}/min_qual_{min_qual}/summary.txt"
    log:
        "logs/truvari/truvari_output.{input_combination}.{min_qual}.log"
    params:
        output_dir = "results/caller_comparison_vaquita_lr/eval/{input_combination}/min_qual_{min_qual}"
    run:
        if (wildcards.input_combination == 'S2') | (wildcards.input_combination == 'L3') | (wildcards.input_combination == 'S2L3'):
            truth_set_gz = config["truth_set_HG002_DEL_renamed_chr"]["gz"],
            truth_set_bed = config["truth_set_HG002_renamed_chr"]["bed"]
        else: # S1, L1, L2, S1L1, S1L2
            truth_set_gz = config["truth_set_HG002_DEL"]["gz"],
            truth_set_bed = config["truth_set_HG002"]["bed"]
        shell("""
        rm -rf {params.output_dir} && truvari bench -b {truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
            --passonly --includebed {truth_set_bed} &>> {log}
        """)
        # => The combinations of S1L3, S2L1, S2L2 are not possible.
        # -f data/reference/hs37d5.fa

rule reformat_truvari_results:
    input:
        "results/caller_comparison_vaquita_lr/eval/{input_combination}/min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_vaquita_lr/eval/{input_combination}/min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.input_combination}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule cat_truvari_results_all:
    input:
        S1   = expand("results/caller_comparison_vaquita_lr/eval/S1/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        # S2 run: (stop after 60h)
        # [2022-4-21 14:51:46] [START] SV CLASSIFICATION
        # [2022-4-21 14:51:46] [START] FIND INVERSIONS
        # [2022-4-21 14:52:15] [END] FIND INVERSIONS (29 seconds.)
        # [2022-4-21 14:52:15] [START] FIND DELETIONS
        # [2022-4-21 14:52:34] [END] FIND DELETIONS (19 seconds.)
        # [2022-4-21 14:52:34] [START] FIND DUPLICATIONS
        # [2022-4-22 2:34:23] [END] FIND DUPLICATIONS (42109 seconds.)
        # [2022-4-22 2:34:23] [START] FIND TRANSLOCATION
        # S2   = expand("results/caller_comparison_vaquita_lr/eval/S2/min_qual_{min_qual}/pr_rec.txt",
        #               min_qual = min_qual_VaquitaLR),
        S1L1 = expand("results/caller_comparison_vaquita_lr/eval/S1L1/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        S1L2 = expand("results/caller_comparison_vaquita_lr/eval/S1L2/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        # S2L3 = expand("results/caller_comparison_vaquita_lr/eval/S2L3/min_qual_{min_qual}/pr_rec.txt",
        #               min_qual = min_qual_VaquitaLR),
        L1   = expand("results/caller_comparison_vaquita_lr/eval/L1/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        L2   = expand("results/caller_comparison_vaquita_lr/eval/L2/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR),
        L3   = expand("results/caller_comparison_vaquita_lr/eval/L3/min_qual_{min_qual}/pr_rec.txt",
                      min_qual = min_qual_VaquitaLR)
    output:
        S1   = temp("results/caller_comparison_vaquita_lr/eval/S1.all_results.txt"),
        # S2   = temp("results/caller_comparison_vaquita_lr/eval/S2.all_results.txt"),
        S1L1 = temp("results/caller_comparison_vaquita_lr/eval/S1L1.all_results.txt"),
        S1L2 = temp("results/caller_comparison_vaquita_lr/eval/S1L2.all_results.txt"),
        # S2L3 = temp("results/caller_comparison_vaquita_lr/eval/S2L3.all_results.txt"),
        L1   = temp("results/caller_comparison_vaquita_lr/eval/L1.all_results.txt"),
        L2   = temp("results/caller_comparison_vaquita_lr/eval/L2.all_results.txt"),
        L3   = temp("results/caller_comparison_vaquita_lr/eval/L3.all_results.txt"),
        all = "results/caller_comparison_vaquita_lr/eval/all_results.txt"
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
