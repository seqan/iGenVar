rule plot_pr_all_results:
    input:
        txt = "results/caller_comparison_short_read/{dataset}/eval/all_results.txt"
    output:
        png = "results/caller_comparison_short_read/{dataset}_results.all.png"
    run:
        if wildcards.dataset == 'Illumina_Paired_End':
            bam_name = config["short_bam_name"]["s1"]
        else: # wildcards.dataset == 'Illumina_Mate_Pair'
            bam_name = config["short_bam_name"]["s2"]
        shell("""
            Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_short_read/workflow/scripts/plot_all_results.R \
                "{bam_name}" {input.txt} Repos/iGenVar/test/benchmark/F1_score.csv {output.png}
        """)
