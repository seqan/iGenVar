rule plot_pr_all_results:
    input:
        txt = "results/caller_comparison_short_read/{dataset}/eval/all_results.txt"
    output:
        png = "results/caller_comparison_short_read/{dataset}_results.all.png"
    run:
        if wildcards.dataset == 'Illumina_Paired_End':
            bam_name = config["short_bam_name"]["s1"]
        elif wildcards.dataset == 'Illumina_Mate_Pair':
            bam_name = config["short_bam_name"]["s2"]
        elif wildcards.dataset == 'hg38_Sim_default':
            bam_name = config["simulated_short_bam_name"]["sim1"]
        elif wildcards.dataset == 'hg38_Sim_InDel':
            bam_name = config["simulated_short_bam_name"]["sim2"]
        elif wildcards.dataset == 'hg38_Sim_noSNP':
            bam_name = config["simulated_short_bam_name"]["sim3"]
        else: # wildcards.dataset == 'hg38_Sim_SNPandSV'
            bam_name = config["simulated_short_bam_name"]["sim4"]
        shell("""
            Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_short_read/workflow/scripts/plot_all_results.R \
                "{bam_name}" {input.txt} Repos/iGenVar/test/benchmark/F1_score.csv {output.png}
        """)
