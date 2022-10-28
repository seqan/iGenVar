min_qual_VaquitaLR = list(range(config["quality_ranges"]["Vaquita-LR"]["from"],
                                config["quality_ranges"]["Vaquita-LR"]["to"],
                                config["quality_ranges"]["Vaquita-LR"]["step"]))

rule simulation_truvari:
    input:
        vcf = "results/caller_comparison_vaquita_lr/{input_combination}/variants.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_vaquita_lr/{input_combination}/variants.min_qual_{min_qual}.vcf.gz.tbi"
    params:
        output_dir = "results/caller_comparison_vaquita_lr/eval/{input_combination}/simulation.min_qual_{min_qual}"
    output:
        summary = "results/caller_comparison_vaquita_lr/eval/{input_combination}/simulation.min_qual_{min_qual}/summary.txt"
    log:
        "logs/truvari/truvari_output.{input_combination}.{min_qual}.log"
    run:
        if (wildcards.input_combination == 'hg38_Sim_default'):
            truth_set_gz = config["truth_set_simulation_default"]["gz"]
        elif (wildcards.input_combination == 'hg38_Sim_InDel'):
            truth_set_gz = config["truth_set_simulation_InDel"]["gz"]
        elif (wildcards.input_combination == 'hg38_Sim_noSNP'):
            truth_set_gz = config["truth_set_simulation_noSNP"]["gz"]
        else: # (wildcards.input_combination == 'hg38_Sim_SNPandSV'):
            truth_set_gz = config["truth_set_simulation_SNPandSV"]["gz"]
        shell("""
            rm -rf {params.output_dir} && truvari bench -b {truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
                --passonly &>> {log}
        """)

rule simulation_reformat_truvari_results:
    input:
        "results/caller_comparison_vaquita_lr/eval/{input_combination}/simulation.min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_vaquita_lr/eval/{input_combination}/simulation.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.input_combination}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule simulation_cat_truvari_results_all:
    input:
        hg38_Sim_default  = expand("results/caller_comparison_vaquita_lr/eval/hg38_Sim_default/simulation.min_qual_{min_qual}/pr_rec.txt",
                                   min_qual = min_qual_VaquitaLR),
        hg38_Sim_InDel    = expand("results/caller_comparison_vaquita_lr/eval/hg38_Sim_InDel/simulation.min_qual_{min_qual}/pr_rec.txt",
                                   min_qual = min_qual_VaquitaLR),
        # hg38_Sim_noSNP    = expand("results/caller_comparison_vaquita_lr/eval/hg38_Sim_noSNP/simulation.min_qual_{min_qual}/pr_rec.txt",
        #                            min_qual = min_qual_VaquitaLR),
        hg38_Sim_SNPandSV = expand("results/caller_comparison_vaquita_lr/eval/hg38_Sim_SNPandSV/simulation.min_qual_{min_qual}/pr_rec.txt",
                                   min_qual = min_qual_VaquitaLR)
    output:
        hg38_Sim_default  = temp("results/caller_comparison_vaquita_lr/eval/hg38_Sim_default.all_results.txt"),
        hg38_Sim_InDel    = temp("results/caller_comparison_vaquita_lr/eval/hg38_Sim_InDel.all_results.txt"),
        # hg38_Sim_noSNP    = temp("results/caller_comparison_vaquita_lr/eval/hg38_Sim_noSNP.all_results.txt"),
        hg38_Sim_SNPandSV = temp("results/caller_comparison_vaquita_lr/eval/hg38_Sim_SNPandSV.all_results.txt"),
        all = "results/caller_comparison_vaquita_lr/eval/simulation.txt"
    threads: 1
    run:
        shell("cat {input.hg38_Sim_default} > {output.hg38_Sim_default}")
        shell("cat {input.hg38_Sim_InDel} > {output.hg38_Sim_InDel}")
        # shell("cat {input.hg38_Sim_noSNP} > {output.hg38_Sim_noSNP}")
        shell("cat {input.hg38_Sim_SNPandSV} > {output.hg38_Sim_SNPandSV}")
        shell("cat {output.hg38_Sim_default} {output.hg38_Sim_InDel} {output.hg38_Sim_SNPandSV} > {output.all}")
        # shell("""
        #     cat {output.hg38_Sim_default} {output.hg38_Sim_InDel} \
        #         {output.hg38_Sim_noSNP} {output.hg38_Sim_SNPandSV} > {output.all}
        # """)
