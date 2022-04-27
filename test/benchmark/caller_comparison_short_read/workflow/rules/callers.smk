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

# Vaquita-LR
# We have only results for S1, S1L1, S1L2
rule copy_vaquita_lr_results:
    output:
        res_S = "results/caller_comparison_short_read/{dataset}/eval/VaquitaLR_S/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
        res_SL1 = "results/caller_comparison_short_read/{dataset}/eval/VaquitaLR_SL1/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
        res_SL2 = "results/caller_comparison_short_read/{dataset}/eval/VaquitaLR_SL2/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
        res_SL3 = "results/caller_comparison_short_read/{dataset}/eval/VaquitaLR_SL3/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    run:
        if wildcards.dataset == 'Illumina_Paired_End':
            shell("""
                cp -r results/caller_comparison_vaquita_lr/eval/S1/DUP_as_INS.min_qual_{wildcards.min_qual} \
                    results/caller_comparison_short_read/{wildcards.dataset}/eval/VaquitaLR_S/
            """)
            shell("""
                cp -r results/caller_comparison_vaquita_lr/eval/S1L1/DUP_as_INS.min_qual_{wildcards.min_qual} \
                    results/caller_comparison_short_read/{wildcards.dataset}/eval/VaquitaLR_SL1/
            """)
            shell("""
                cp -r results/caller_comparison_vaquita_lr/eval/S1L2/DUP_as_INS.min_qual_{wildcards.min_qual} \
                    results/caller_comparison_short_read/{wildcards.dataset}/eval/VaquitaLR_SL2/
            """)
            # As Vaquita LR S1L3 does not exist, we create empty files.
            shell("""
                # cp -r results/caller_comparison_vaquita_lr/eval/S1L3/DUP_as_INS.min_qual_{wildcards.min_qual} \
                #     results/caller_comparison_short_read/{wildcards.dataset}/eval/VaquitaLR_SL3/
                echo 'VaquitaLR_SL3\t{wildcards.min_qual}\tprecision\t0\nVaquitaLR_SL3\t{wildcards.min_qual}\trecall\t0' > {output.res_SL3}
            """)
            shell("""
                sed -i 's/S1/VaquitaLR_S/g' {output.res_S}
                sed -i 's/S1L1/VaquitaLR_SL1/g' {output.res_SL1}
                sed -i 's/S1L2/VaquitaLR_SL2/g' {output.res_SL2}
                sed -i 's/S1L3/VaquitaLR_SL3/g' {output.res_SL3}
            """)
        else: # wildcards.dataset == 'Illumina_Mate_Pair'
            # As Vaquita LR S2 and all combinations not exist, we create empty files.
            shell("""
                # cp -r results/caller_comparison_vaquita_lr/eval/S2/DUP_as_INS.min_qual_{wildcards.min_qual} \
                #     results/caller_comparison_short_read/{wildcards.dataset}/eval/VaquitaLR_S/
                echo 'VaquitaLR_S\t{wildcards.min_qual}\tprecision\t0\nVaquitaLR_S\t{wildcards.min_qual}\trecall\t0' > {output.res_S}
            """)
            shell("""
                # cp -r results/caller_comparison_vaquita_lr/eval/S2L1/DUP_as_INS.min_qual_{wildcards.min_qual} \
                #     results/caller_comparison_short_read/{wildcards.dataset}/eval/VaquitaLR_SL1/
                echo 'VaquitaLR_SL1\t{wildcards.min_qual}\tprecision\t0\nVaquitaLR_SL1\t{wildcards.min_qual}\trecall\t0' > {output.res_SL1}
            """)
            shell("""
                # cp -r results/caller_comparison_vaquita_lr/eval/S2L2/DUP_as_INS.min_qual_{wildcards.min_qual} \
                #     results/caller_comparison_short_read/{wildcards.dataset}/eval/VaquitaLR_SL2/
                echo 'VaquitaLR_SL2\t{wildcards.min_qual}\tprecision\t0\nVaquitaLR_SL2\t{wildcards.min_qual}\trecall\t0' > {output.res_SL2}
            """)
            shell("""
                # cp -r results/caller_comparison_vaquita_lr/eval/S2L3/DUP_as_INS.min_qual_{wildcards.min_qual} \
                #     results/caller_comparison_short_read/{wildcards.dataset}/eval/VaquitaLR_SL3/
                echo 'VaquitaLR_SL3\t{wildcards.min_qual}\tprecision\t0\nVaquitaLR_SL3\t{wildcards.min_qual}\trecall\t0' > {output.res_SL3}
            """)
            shell("""
                sed -i 's/S2/VaquitaLR_S/g' {output.res_S}
                sed -i 's/S2L1/VaquitaLR_SL1/g' {output.res_SL1}
                sed -i 's/S2L2/VaquitaLR_SL2/g' {output.res_SL2}
                sed -i 's/S2L3/VaquitaLR_SL3/g' {output.res_SL3}
            """)
