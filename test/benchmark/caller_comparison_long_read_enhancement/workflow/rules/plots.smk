rule plot_pr_all_results:
    input:
        a = "results/caller_comparison_long_read_enhancement/{dataset}/eval/all_results.txt",
        b = "results/caller_comparison_long_read_enhancement/{dataset}/eval/DUP_as_INS.all_results.txt"
    output:
        a = "results/caller_comparison_long_read_enhancement/{dataset}/eval/results.all.png",
        b = "results/caller_comparison_long_read_enhancement/{dataset}/eval/results.DUP_as_INS.all.png"
    run:
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_long_read_enhancement/workflow/scripts/plot_all_results.R \
            {input.a} Repos/iGenVar/test/benchmark/F1_score.csv {output.a}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_long_read_enhancement/workflow/scripts/plot_all_results.R \
            {input.b} Repos/iGenVar/test/benchmark/F1_score.csv {output.b}
        """)
