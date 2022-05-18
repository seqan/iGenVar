rule reformat_truvari_results_2:
    input:
        "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/min_qual_{min_qual}/summary.txt"
    output:
        pr = "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/min_qual_{min_qual}/pr.txt",
        rec = "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/min_qual_{min_qual}/rec.txt"
    run:
        shell("""
            cat {input} | grep '\<TP-base\>\|\<FP\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr -d "-" | tr -d base \
                | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.long_read_enhancement}\", \"{wildcards.min_qual}\", $1, $2 }}' \
                > {output.pr}
        """)
        shell("""
            cat {input} | grep '\<TP-base\>\|\<FN\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr -d "-" | tr -d base \
                | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.long_read_enhancement}\", \"{wildcards.min_qual}\", $1, $2 }}' \
                > {output.rec}
        """)

rule cat_truvari_pr_rec_results_all:
    input:
        pr_S      = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/no_enhancement/min_qual_{min_qual}/pr.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        pr_SL2x1  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x1/min_qual_{min_qual}/pr.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        pr_SL2x2  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x2/min_qual_{min_qual}/pr.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        pr_SL2x3  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x3/min_qual_{min_qual}/pr.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        pr_L2x1   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x1/min_qual_{min_qual}/pr.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        pr_L2x2   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x2/min_qual_{min_qual}/pr.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        pr_L2x3   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x3/min_qual_{min_qual}/pr.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        pr_L2     = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2/min_qual_{min_qual}/pr.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        rec_S     = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/no_enhancement/min_qual_{min_qual}/rec.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        rec_SL2x1 = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x1/min_qual_{min_qual}/rec.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        rec_SL2x2 = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x2/min_qual_{min_qual}/rec.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        rec_SL2x3 = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x3/min_qual_{min_qual}/rec.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        rec_L2x1  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x1/min_qual_{min_qual}/rec.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        rec_L2x2  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x2/min_qual_{min_qual}/rec.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        rec_L2x3  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x3/min_qual_{min_qual}/rec.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
        rec_L2    = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2/min_qual_{min_qual}/rec.txt",
                       min_qual=list(range(config["quality_ranges"]["iGenVar"]["from"],
                                           config["quality_ranges"]["iGenVar"]["to"],
                                           config["quality_ranges"]["iGenVar"]["step"]))),
    output:
        pr_SL2x1 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/pr_SL2x1.all_results.txt"),
        pr_S     = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/pr_S.all_results.txt"),
        pr_SL2x2 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/pr_SL2x2.all_results.txt"),
        pr_SL2x3 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/pr_SL2x3.all_results.txt"),
        pr_L2x1 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/pr_L2x1.all_results.txt"),
        pr_L2x2 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/pr_L2x2.all_results.txt"),
        pr_L2x3 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/pr_L2x3.all_results.txt"),
        pr_L2    = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/pr_L2.all_results.txt"),
        rec_S     = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/rec_S.all_results.txt"),
        rec_SL2x1 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/rec_SL2x1.all_results.txt"),
        rec_SL2x2 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/rec_SL2x2.all_results.txt"),
        rec_SL2x3 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/rec_SL2x3.all_results.txt"),
        rec_L2x1 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/rec_L2x1.all_results.txt"),
        rec_L2x2 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/rec_L2x2.all_results.txt"),
        rec_L2x3 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/rec_L2x3.all_results.txt"),
        rec_L2    = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/rec_L2.all_results.txt"),
        pr_all = "results/caller_comparison_long_read_enhancement/{dataset}/eval/pr_all_results.txt",
        rec_all = "results/caller_comparison_long_read_enhancement/{dataset}/eval/rec_all_results.txt"
    run:
        shell("cat {input.pr_S} > {output.pr_S}")
        shell("cat {input.pr_SL2x1} > {output.pr_SL2x1}")
        shell("cat {input.pr_SL2x2} > {output.pr_SL2x2}")
        shell("cat {input.pr_SL2x3} > {output.pr_SL2x3}")
        shell("cat {input.pr_L2x1} > {output.pr_L2x1}")
        shell("cat {input.pr_L2x2} > {output.pr_L2x2}")
        shell("cat {input.pr_L2x3} > {output.pr_L2x3}")
        shell("cat {input.pr_L2} > {output.pr_L2}")
        shell("""
            cat {output.pr_S} {output.pr_SL2x1} {output.pr_SL2x2} {output.pr_SL2x3} \
                {output.pr_L2x1} {output.pr_L2x2} {output.pr_L2x3} {output.pr_L2} > {output.pr_all}
        """)
        shell("cat {input.rec_S} > {output.rec_S}")
        shell("cat {input.rec_SL2x1} > {output.rec_SL2x1}")
        shell("cat {input.rec_SL2x2} > {output.rec_SL2x2}")
        shell("cat {input.rec_SL2x3} > {output.rec_SL2x3}")
        shell("cat {input.rec_L2x1} > {output.rec_L2x1}")
        shell("cat {input.rec_L2x2} > {output.rec_L2x2}")
        shell("cat {input.rec_L2x3} > {output.rec_L2x3}")
        shell("cat {input.rec_L2} > {output.rec_L2}")
        shell("""
            cat {output.rec_S} {output.rec_SL2x1} {output.rec_SL2x2} {output.rec_SL2x3} \
                {output.rec_L2x1} {output.rec_L2x2} {output.rec_L2x3} {output.rec_L2} > {output.rec_all}
        """)
