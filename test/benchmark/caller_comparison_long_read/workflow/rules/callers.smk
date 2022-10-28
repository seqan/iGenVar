sample = config["parameters"]["sample"],
min_var_length = config["parameters"]["min_var_length"],
max_var_length = config["parameters"]["max_var_length"]

# iGenVar
rule copy_igenvar_results:
    output:
        res_L = "results/caller_comparison_long_read/{dataset}/eval/iGenVar_L/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
        res_SL = "results/caller_comparison_long_read/{dataset}/eval/iGenVar_SL/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    run:
        if wildcards.dataset == 'MtSinai_PacBio':
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/L1/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_long_read/MtSinai_PacBio/eval/iGenVar_L/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S1L1/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_long_read/MtSinai_PacBio/eval/iGenVar_SL/
            """)
            shell("""
                sed -i 's/L1/iGenVar_L/g' {output.res_L}
                sed -i 's/S1L1/iGenVar_SL/g' {output.res_SL}
            """)
        elif wildcards.dataset == 'PacBio_CCS':
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/L2/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_long_read/PacBio_CCS/eval/iGenVar_L/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S1L2/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_long_read/PacBio_CCS/eval/iGenVar_SL/
            """)
            shell("""
                sed -i 's/L2/iGenVar_L/g' {output.res_L}
                sed -i 's/S1L2/iGenVar_SL/g' {output.res_SL}
            """)
        else: # wildcards.dataset == '10X_Genomics'
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/L3/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_long_read/10X_Genomics/eval/iGenVar_L/
            """)
            shell("""
            cp -r results/caller_comparison_iGenVar_only/eval/S1L3/no_DUP_and_INV.min_qual_{wildcards.min_qual} \
                results/caller_comparison_long_read/10X_Genomics/eval/iGenVar_SL/
            """)
            shell("""
                sed -i 's/L3/iGenVar_L/g' {output.res_L}
                sed -i 's/S1L3/iGenVar_SL/g' {output.res_SL}
            """)

# SVIM
rule run_svim:
    output:
        vcf = "results/caller_comparison_long_read/{dataset}/SVIM/unsorted/variants.vcf"
    log:
        "logs/caller_comparison_long_read/SVIM/svim_output.{dataset}.log"
    params:
        working_dir = "results/caller_comparison_long_read/{dataset}/SVIM/unsorted/"
    threads: 1
    run:
        # if wildcards.dataset == 'MtSinai_PacBio':
        #     long_bam = config["long_read_bam"]["l1"],
        #     genome = config["reference_fa"]["MtSinai_PacBio"]
        if wildcards.dataset == 'PacBio_CCS':
            long_bam = config["long_read_bam"]["l2"],
            genome = config["reference_fa"]["PacBio_CCS"]
        else: # wildcards.dataset == '10X_Genomics'
            long_bam = config["long_read_bam"]["l3"],
            genome = config["reference_fa"]["10X_Genomics"]
        shell("""
            /usr/bin/time -v svim alignment --sample {sample} --partition_max_distance 1000 --cluster_max_distance 0.5 \
                --min_sv_size {min_var_length} --segment_gap_tolerance 20 --segment_overlap_tolerance 20 \
                --read_names --max_sv_size {max_var_length} {params.working_dir} {long_bam} {genome} &>> {log}
        """)
        # Defaults:
        # --position_distance_normalizer 900 --edit_distance_normalizer 1.0

# SVIM does not output sorted vcf files, thus we have to sort them via picard.
rule picard:
    input:
        vcf = "results/caller_comparison_long_read/{dataset}/SVIM/unsorted/variants.vcf"
    output:
        vcf = "results/caller_comparison_long_read/{dataset}/SVIM/variants.vcf"
    log:
        "logs/caller_comparison_long_read/SVIM/picard_output.{dataset}.log"
    conda:
        "../../../envs/simulation.yaml"
    shell:
        "picard SortVcf -I {input.vcf} -O {output.vcf} -Xms1g -Xmx100g --TMP_DIR tmp/picard/ &>> {log}"
        # The Xms and Xmx sets the java memory for avoiding "java.lang.OutOfMemoryError: GC overhead limit exceeded"

# As SVIM is not working with the dataset MtSinai_PacBio we create an empty output file.
rule create_empty_SVIM:
    output:
        txt = "results/caller_comparison_long_read/MtSinai_PacBio/eval/SVIM/min_qual_{min_qual}/pr_rec.txt"
    shell:
        "echo 'SVIM\t{wildcards.min_qual}\tprecision\t0\nSVIM\t{wildcards.min_qual}\trecall\t0' > {output.txt}"

# Vaquita-LR
# We have only results for L1, L2, L3, S1L1, S1L2
rule copy_vaquita_lr_results:
    output:
        res_L = "results/caller_comparison_long_read/{dataset}/eval/Vaquita_lr_L/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
        res_SL = "results/caller_comparison_long_read/{dataset}/eval/Vaquita_lr_SL/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    run:
        if wildcards.dataset == 'MtSinai_PacBio':
            shell("""
                cp -r results/caller_comparison_vaquita_lr/eval/L1/DUP_as_INS.min_qual_{wildcards.min_qual} \
                    results/caller_comparison_long_read/{wildcards.dataset}/eval/Vaquita_lr_L/
            """)
            shell("""
                cp -r results/caller_comparison_vaquita_lr/eval/S1L1/DUP_as_INS.min_qual_{wildcards.min_qual} \
                    results/caller_comparison_long_read/{wildcards.dataset}/eval/Vaquita_lr_SL/
            """)
            shell("""
                sed -i 's/L1/Vaquita_lr_L/g' {output.res_L}
                sed -i 's/S1L1/Vaquita_lr_SL/g' {output.res_SL}
            """)
        elif wildcards.dataset == 'PacBio_CCS':
            shell("""
                cp -r results/caller_comparison_vaquita_lr/eval/L2/DUP_as_INS.min_qual_{wildcards.min_qual} \
                    results/caller_comparison_long_read/{wildcards.dataset}/eval/Vaquita_lr_L/
            """)
            shell("""
                cp -r results/caller_comparison_vaquita_lr/eval/S1L2/DUP_as_INS.min_qual_{wildcards.min_qual} \
                    results/caller_comparison_long_read/{wildcards.dataset}/eval/Vaquita_lr_SL/
            """)
            shell("""
                sed -i 's/L2/Vaquita_lr_L/g' {output.res_L}
                sed -i 's/S1L2/Vaquita_lr_SL/g' {output.res_SL}
            """)
        else: # wildcards.dataset == '10X_Genomics'
            shell("""
            cp -r results/caller_comparison_vaquita_lr/eval/L3/DUP_as_INS.min_qual_{wildcards.min_qual} \
                results/caller_comparison_long_read/{wildcards.dataset}/eval/Vaquita_lr_L/
            """)
            # As Vaquita LR S1L3 does not exist, we create empty files.
            shell("""
                echo 'Vaquita_lr_SL\t{wildcards.min_qual}\tprecision\t0\nVaquita_lr_SL\t{wildcards.min_qual}\trecall\t0' > {output.res_SL}
            """)
            shell("""
                sed -i 's/L3/Vaquita_lr_L/g' {output.res_L}
            """)

# SNIFFLES (we have to loop over min_support, because sniffles does not write a quality score into the vcf)
rule run_sniffles:
    output:
        "results/caller_comparison_long_read/{dataset}/Sniffles/raw_variants_{min_qual}.vcf"
    log:
        "logs/caller_comparison_long_read/sniffles/sniffles_output.{dataset}.{min_qual}.log"
    threads: 8
    run:
        if wildcards.dataset == 'MtSinai_PacBio':
            long_bam = config["long_read_md_bam"]["l1"]
        elif wildcards.dataset == 'PacBio_CCS':
            long_bam = config["long_read_md_bam"]["l2"]
        else: # wildcards.dataset == '10X_Genomics'
            long_bam = config["long_read_md_bam"]["l3"]
        shell("""
            /usr/bin/time -v sniffles --mapped_reads {long_bam} --vcf {output} \
                --min_support {wildcards.min_qual} --min_length {min_var_length} --threads {threads} --genotype &>> {log}
        """)

# see https://github.com/spiralgenetics/truvari/issues/43
# see https://github.com/fritzsedlazeck/Sniffles/issues/209 (Fixed, but there just hasn't been a release since the fix.)
rule fix_sniffles:
    input:
        "results/caller_comparison_long_read/{dataset}/Sniffles/raw_variants_{min_qual}.vcf"
    output:
        "results/caller_comparison_long_read/{dataset}/Sniffles/variants.unsorted.min_qual_{min_qual}.vcf"
    run:
        shell("sed 's/##INFO=<ID=SUPTYPE,Number=A/##INFO=<ID=SUPTYPE,Number=./' {input} > {output}")
        shell("sed -i '4i##FILTER=<ID=STRANDBIAS,Description=\"Strand is biased.\">' {output}")

# Split to SV classes
rule fix_sniffles_2:
    input:
        "results/caller_comparison_long_read/{dataset}/Sniffles/variants.unsorted.min_qual_{min_qual}.vcf"
    output:
        "results/caller_comparison_long_read/{dataset}/Sniffles/variants.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools sort {input} > {output}"

#PBSV
rule run_pbsv_dicsover:
    output:
        svsig_gz = "results/caller_comparison_long_read/{dataset}/pbsv/signatures.svsig.gz"
        # svsig_gz = dynamic("results/caller_comparison_long_read/{dataset}/pbsv/signatures/{region}.svsig.gz")
    log:
        log = "logs/caller_comparison_long_read/pbsv/pbsv_discover.{dataset}.log"
    threads: 8
    run:
        if wildcards.dataset == 'MtSinai_PacBio':
            long_bam = config["long_read_bam"]["l1"]
        elif wildcards.dataset == 'PacBio_CCS':
            long_bam = config["long_read_bam"]["l2"]
        else: # wildcards.dataset == '10X_Genomics'
            long_bam = config["long_read_bam"]["l3"]
        shell("/usr/bin/time -v pbsv discover {long_bam} {output.svsig_gz} &>> {log}")
        # shell("""
        #     for i in $(samtools view -H {long_bam} | grep '^@SQ' | cut -f2 | cut -d':' -f2); do
        #         /usr/bin/time -v pbsv discover --region $i {long_bam} results/caller_comparison_long_read/{wildcards.dataset}/pbsv/signatures/$i.svsig.gz &>> {output.log} &
        #     done
        #     wait
        # """)

# pbsv call is not using more than 200% CPU capacity (2 threads)
rule run_pbsv_call:
    input:
        svsig_gz = "results/caller_comparison_long_read/{dataset}/pbsv/signatures.svsig.gz"
    output:
        vcf = "results/caller_comparison_long_read/{dataset}/pbsv/variants.min_qual_{min_qual}.vcf"
    log:
        "logs/caller_comparison_long_read/pbsv/pbsv_call.{dataset}.{min_qual}.log"
    threads: 2
    run:
        if wildcards.dataset == 'MtSinai_PacBio':
            genome = config["reference_fa"]["MtSinai_PacBio"]
        elif wildcards.dataset == 'PacBio_CCS':
            genome = config["reference_fa"]["PacBio_CCS"]
        else: # wildcards.dataset == '10X_Genomics'
            genome = config["reference_fa"]["10X_Genomics"]
        shell("""
            /usr/bin/time -v pbsv call --types DEL,INS,DUP --min-sv-length {min_var_length} --max-ins-length 100K \
                --call-min-reads-all-samples {wildcards.min_qual} --call-min-reads-one-sample {wildcards.min_qual} \
                --call-min-reads-per-strand-all-samples 0 --call-min-bnd-reads-all-samples 0 --call-min-read-perc-one-sample 0 \
                --num-threads {threads} {genome} {input.svsig_gz} {output.vcf} &>> {log}
        """)

# pbsv call is not using more than 200% CPU capacity (2 threads)
rule run_pbsv_call_without_DUP:
    input:
        svsig_gz = "results/caller_comparison_long_read/{dataset}/pbsv/signatures.svsig.gz"
    output:
        vcf = "results/caller_comparison_long_read/{dataset}/pbsv_without_DUP/variants.min_qual_{min_qual}.vcf"
    log:
        "logs/caller_comparison_long_read/pbsv/pbsv_call_without_DUP.{dataset}.{min_qual}.log"
    threads: 2
    run:
        if wildcards.dataset == 'MtSinai_PacBio':
            genome = config["reference_fa"]["MtSinai_PacBio"]
        elif wildcards.dataset == 'PacBio_CCS':
            genome = config["reference_fa"]["PacBio_CCS"]
        else: # wildcards.dataset == '10X_Genomics'
            genome = config["reference_fa"]["10X_Genomics"]
        shell("""
            /usr/bin/time -v pbsv call --types DEL,INS --min-sv-length {min_var_length} --max-ins-length 100K \
                --call-min-reads-all-samples {wildcards.min_qual} --call-min-reads-one-sample {wildcards.min_qual} \
                --call-min-reads-per-strand-all-samples 0 --call-min-bnd-reads-all-samples 0 --call-min-read-perc-one-sample 0 \
                --num-threads {threads} {genome} {input.svsig_gz} {output.vcf} &>> {log}
        """)
