rule plot_pr_all_results:
    input:
        "results/caller_comparison_iGenVar_only/eval/all_results.txt"
    output:
        "results/caller_comparison_iGenVar_only/eval/results.all.png"
    log:
        "logs/caller_comparison_iGenVar_only/rplot.all.log"
    shell:
        """
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            {input} Repos/iGenVar/test/benchmark/F1_score.csv {output} > {log}
        """
