rule plot_pr_all_results:
    input:
        "results/caller_comparison/eval/all_results.txt"
    output:
        "results/caller_comparison/eval/results.all.png"
    log:
        "logs/rplot.all.log"
    shell:
        """
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison/workflow/scripts/plot_all_results.R \
        {input} Repos/iGenVar/test/benchmark/F1_score.csv {output} > {log}
        """
