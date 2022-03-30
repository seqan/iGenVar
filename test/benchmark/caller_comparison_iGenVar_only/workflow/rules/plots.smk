rule plot_pr_all_results:
    input:
        a = "results/caller_comparison_iGenVar_only/eval/all_results.txt",
        b = "results/caller_comparison_iGenVar_only/eval/DUP_as_INS.all_results.txt"
    output:
        a = "results/caller_comparison_iGenVar_only/eval/results.all.png",
        b = "results/caller_comparison_iGenVar_only/eval/results.DUP_as_INS.all.png"
    log:
        a = "logs/caller_comparison_iGenVar_only/rplot.all.log",
        b = "logs/caller_comparison_iGenVar_only/rplot.DUP_as_INS.all.log"
    run:
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            {input.a} Repos/iGenVar/test/benchmark/F1_score.csv {output.a} > {log.a}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            {input.b} Repos/iGenVar/test/benchmark/F1_score.csv {output.b} > {log.b}
        """)
