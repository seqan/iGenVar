rule plot_pr_all_results:
    input:
        "results/caller_comparison_long_read/{dataset}/eval/all_results.txt"
    output:
        "results/caller_comparison_long_read/{dataset}_results.all.png"
    run:
        if wildcards.dataset == 'MtSinai_PacBio':
            bam_name = config["long_bam_name"]["l1"]
        elif wildcards.dataset == 'PacBio_CCS':
            bam_name = config["long_bam_name"]["l2"]
        else: # wildcards.dataset == '10X_Genomics'
            bam_name = config["long_bam_name"]["l3"]
        shell("""
            Rscript --vanilla Repos/iGenVar/test/benchmark/caller_comparison_long_read/workflow/scripts/plot_all_results.R \
                "{bam_name}" {input} Repos/iGenVar/test/benchmark/F1_score.csv {output}
        """)
