sample = config["parameters"]["sample"],
min_var_length = config["parameters"]["min_var_length"],
max_var_length = config["parameters"]["max_var_length"]

# iGenVar
rule copy_igenvar_results:
    output:
        res_S = "results/caller_comparison_short_read/{dataset}/eval/iGenVar_S/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
        res_SL1 = "results/caller_comparison_short_read/{dataset}/eval/iGenVar_SL1/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
        res_SL2 = "results/caller_comparison_short_read/{dataset}/eval/iGenVar_SL2/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
        res_SL3 = "results/caller_comparison_short_read/{dataset}/eval/iGenVar_SL3/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    run:
        if wildcards.dataset == 'Illumina_Paired_End':
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S1/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_short_read/Illumina_Paired_End/eval/iGenVar_S/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S1L1/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_short_read/Illumina_Paired_End/eval/iGenVar_SL1/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S1L2/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_short_read/Illumina_Paired_End/eval/iGenVar_SL2/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S1L3/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_short_read/Illumina_Paired_End/eval/iGenVar_SL3/
            """)
            shell("""
                sed -i 's/S1/iGenVar_S/g' {output.res_S}
                sed -i 's/S1L1/iGenVar_SL1/g' {output.res_SL1}
                sed -i 's/S1L2/iGenVar_SL2/g' {output.res_SL2}
                sed -i 's/S1L3/iGenVar_SL3/g' {output.res_SL3}
            """)
        else: # wildcards.dataset == 'Illumina_Mate_Pair'
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S2/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_short_read/Illumina_Mate_Pair/eval/iGenVar_S/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S2L1/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_short_read/Illumina_Mate_Pair/eval/iGenVar_SL1/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S2L2/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_short_read/Illumina_Mate_Pair/eval/iGenVar_SL2/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S2L3/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_short_read/Illumina_Mate_Pair/eval/iGenVar_SL3/
            """)
            shell("""
                sed -i 's/S2/iGenVar_S/g' {output.res_S}
                sed -i 's/S2L1/iGenVar_SL1/g' {output.res_SL1}
                sed -i 's/S2L2/iGenVar_SL2/g' {output.res_SL2}
                sed -i 's/S2L3/iGenVar_SL3/g' {output.res_SL3}
            """)
