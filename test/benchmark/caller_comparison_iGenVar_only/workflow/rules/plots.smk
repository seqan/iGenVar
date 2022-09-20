rule plot_pr_all_results:
    input:
        a = "results/caller_comparison_iGenVar_only/eval/all_results.txt",
        b = "results/caller_comparison_iGenVar_only/eval/DEL_results.txt",
        c = "results/caller_comparison_iGenVar_only/eval/INS_results.txt",
        d = "results/caller_comparison_iGenVar_only/eval/no_DUP_and_INV_results.txt"
    output:
        a = "results/caller_comparison_iGenVar_only/results.all.png",
        b = "results/caller_comparison_iGenVar_only/results.DEL.png",
        c = "results/caller_comparison_iGenVar_only/results.INS.png",
        d = "results/caller_comparison_iGenVar_only/results.no_DUP_and_INV.png"
    run:
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            "iGenVar 0.0.3" {input.a} Repos/iGenVar/test/benchmark/F1_score.csv {output.a}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            "iGenVar 0.0.3 - Deletion only" {input.b} Repos/iGenVar/test/benchmark/F1_score.csv {output.b}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            "iGenVar 0.0.3 - Insertion only" {input.c} Repos/iGenVar/test/benchmark/F1_score.csv {output.c}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_iGenVar_only/workflow/scripts/plot_all_results.R \
            "iGenVar 0.0.3 - Duplications and Inversions as Insertions" {input.d} Repos/iGenVar/test/benchmark/F1_score.csv {output.d}
        """)
