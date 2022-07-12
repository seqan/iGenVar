rule plot_pr_all_results:
    input:
        a = "results/caller_comparison_iGenVar_only/eval/all_results.txt",
        b = "results/caller_comparison_iGenVar_only/eval/no_DUP_and_INV.all_results.txt"
    output:
        a = "results/caller_comparison_iGenVar_only/eval/results.all.png",
        b = "results/caller_comparison_iGenVar_only/eval/results.no_DUP_and_INV.all.png"
    log:
        a = "logs/caller_comparison_iGenVar_only/rplot.all.log",
        b = "logs/caller_comparison_iGenVar_only/rplot.no_DUP_and_INV.all.log"
    run:
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            {input.a} Repos/iGenVar/test/benchmark/F1_score.csv {output.a} > {log.a}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            {input.b} Repos/iGenVar/test/benchmark/F1_score.csv {output.b} > {log.b}
        """)
