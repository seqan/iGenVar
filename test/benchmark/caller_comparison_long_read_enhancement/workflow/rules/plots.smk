rule plot_all_results:
    input:
        a = "results/caller_comparison_long_read_enhancement/{dataset}/eval/all_results.txt",
        b = "results/caller_comparison_long_read_enhancement/{dataset}/eval/DUP_as_INS.all_results.txt"
    output:
        a = "results/caller_comparison_long_read_enhancement/{dataset}.results.all.png",
        b = "results/caller_comparison_long_read_enhancement/{dataset}.results.DUP_as_INS.all.png"
    run:
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_long_read_enhancement/workflow/scripts/plot_all_results.R \
            {input.a} Repos/iGenVar/test/benchmark/F1_score.csv {output.a}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_long_read_enhancement/workflow/scripts/plot_all_results.R \
            {input.b} Repos/iGenVar/test/benchmark/F1_score.csv {output.b}
        """)

rule plot_pr_rec_all_results:
    input:
        pr = "results/caller_comparison_long_read_enhancement/{dataset}/eval/pr_all_results.txt",
        rec = "results/caller_comparison_long_read_enhancement/{dataset}/eval/rec_all_results.txt"
    output:
        pr = "results/caller_comparison_long_read_enhancement/{dataset}.results.pr_all.png",
        rec = "results/caller_comparison_long_read_enhancement/{dataset}.results.rec_all.png"
    run:
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_long_read_enhancement/workflow/scripts/plot_pr_all_results.R \
            {input.pr} Repos/iGenVar/test/benchmark/Precision_or_Recall_score.csv {output.pr}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_long_read_enhancement/workflow/scripts/plot_rec_all_results.R \
            {input.rec} Repos/iGenVar/test/benchmark/Precision_or_Recall_score.csv {output.rec}
        """)
