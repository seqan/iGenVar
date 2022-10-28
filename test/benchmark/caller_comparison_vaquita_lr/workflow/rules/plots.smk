rule plot_pr_all_results:
    input:
        a = "results/caller_comparison_vaquita_lr/eval/all_results.txt",
        b = "results/caller_comparison_vaquita_lr/eval/DUP_as_INS.all_results.txt",
        c = "results/caller_comparison_vaquita_lr/eval/simulation.txt"
    output:
        a = "results/caller_comparison_vaquita_lr/results.DEL_only.png",
        b = "results/caller_comparison_vaquita_lr/results.DEL_and_DUP_as_INS.png",
        c = "results/caller_comparison_vaquita_lr/results.simulation.png"
    run:
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_vaquita_lr/workflow/scripts/plot_all_results.R \
            "Vaquita LR - Deletions only" {input.a} Repos/iGenVar/test/benchmark/F1_score.csv {output.a}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_vaquita_lr/workflow/scripts/plot_all_results.R \
            "Vaquita LR - Deletions + Duplications as Insertions" {input.b} Repos/iGenVar/test/benchmark/F1_score.csv {output.b}
        """)
        shell("""
        Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_vaquita_lr/workflow/scripts/simulation.plot_all_results.R \
            {input.c} Repos/iGenVar/test/benchmark/F1_score.csv {output.c}
        """)
