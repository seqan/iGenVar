min_qual_iGenVar   = list(range(config["quality_ranges"]["iGenVar"]["from"],
                                config["quality_ranges"]["iGenVar"]["to"],
                                config["quality_ranges"]["iGenVar"]["step"]))
min_qual_VaquitaLR = list(range(config["quality_ranges"]["Vaquita-LR"]["from"],
                                config["quality_ranges"]["Vaquita-LR"]["to"],
                                config["quality_ranges"]["Vaquita-LR"]["step"]))

rule filter_vaquita_vcf:
    input:
        vcf = "results/caller_comparison_short_read/{dataset}/Vaquita/variants.vcf"
    output:
        vcf = "results/caller_comparison_short_read/{dataset}/Vaquita/variants.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'INFO/SC>={wildcards.min_qual}' {input.vcf} > {output.vcf}"

rule filter_bcf:
    input:
        bcf = "results/caller_comparison_short_read/{dataset}/Delly2/variants.bcf"
    output:
        vcf ="results/caller_comparison_short_read/{dataset}/Delly2/variants.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools convert {input.bcf} | bcftools view -i 'QUAL>={wildcards.min_qual}' > {output.vcf}"

rule filter_vcf:
    input:
        vcf = "results/caller_comparison_short_read/{dataset}/{caller,GRIDSS|TIDDIT}/variants.vcf"
    output:
        vcf = "results/caller_comparison_short_read/{dataset}/{caller,GRIDSS|TIDDIT}/variants.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'QUAL>={wildcards.min_qual}' {input.vcf} > {output.vcf}"

rule bgzip:
    input:
        "{name}.vcf"
    output:
        "{name}.vcf.gz"
    shell:
        "bgzip -c {input} > {output}"

rule tabix:
    input:
        "{name}.vcf.gz"
    output:
        "{name}.vcf.gz.tbi"
    shell:
        "tabix -p vcf {input}"

rule truvari:
    input:
        vcf = "results/caller_comparison_short_read/{dataset}/{caller}/variants.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_short_read/{dataset}/{caller}/variants.min_qual_{min_qual}.vcf.gz.tbi"
    params:
        output_dir = "results/caller_comparison_short_read/{dataset}/eval/{caller}/min_qual_{min_qual}"
    output:
        summary = "results/caller_comparison_short_read/{dataset}/eval/{caller}/min_qual_{min_qual}/summary.txt"
    log:
        "logs/caller_comparison_short_read/truvari/truvari_output.{dataset}.{caller}.{min_qual}.log"
    run:
        if (wildcards.dataset == 'hg38_Sim_default'):
            truth_set_gz = config["truth_set_simulation_default"]["gz"]
        elif (wildcards.dataset == 'hg38_Sim_InDel'):
            truth_set_gz = config["truth_set_simulation_InDel"]["gz"]
        elif (wildcards.dataset == 'hg38_Sim_noSNP'):
            truth_set_gz = config["truth_set_simulation_noSNP"]["gz"]
        elif (wildcards.dataset == 'hg38_Sim_SNPandSV'):
            truth_set_gz = config["truth_set_simulation_SNPandSV"]["gz"]
        elif (wildcards.caller == 'GRIDSS'):
            truth_set_gz = config["truth_set_HG002_renamed_chr"]["gz"],
            truth_set_bed = config["truth_set_HG002_renamed_chr"]["bed"]
        elif wildcards.dataset == 'Illumina_Paired_End':
            truth_set_gz = config["truth_set_HG002"]["gz"],
            truth_set_bed = config["truth_set_HG002"]["bed"]
        else: # wildcards.dataset == 'Illumina_Mate_Pair'
            truth_set_gz = config["truth_set_HG002_renamed_chr"]["gz"],
            truth_set_bed = config["truth_set_HG002_renamed_chr"]["bed"]
        if (wildcards.dataset == 'Illumina_Paired_End') | (wildcards.dataset == 'Illumina_Mate_Pair'):
            shell("""
                rm -rf {params.output_dir} && truvari bench -b {truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
                    --passonly --includebed {truth_set_bed} &>> {log}
            """)
        else: # wildcards.dataset == 'hg38_Sim_*'
            shell("""
                rm -rf {params.output_dir} && truvari bench -b {truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
                    --passonly &>> {log}
            """)

rule reformat_truvari_results:
    input:
        txt = "results/caller_comparison_short_read/{dataset}/eval/{caller}/min_qual_{min_qual}/summary.txt"
    output:
        txt = "results/caller_comparison_short_read/{dataset}/eval/{caller}/min_qual_{min_qual}/pr_rec.txt"
    run:
        shell("""
            cat {input.txt} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
                | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.min_qual}\", $1, $2 }}' \
                > {output.txt}
        """)

rule cat_truvari_results_all:
    input:
        iGenVar_S     = expand("results/caller_comparison_short_read/{{dataset}}/eval/iGenVar_S/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
                               min_qual = min_qual_iGenVar),
        iGenVar_SL1   = expand("results/caller_comparison_short_read/{{dataset}}/eval/iGenVar_SL1/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
                               min_qual = min_qual_iGenVar),
        iGenVar_SL2   = expand("results/caller_comparison_short_read/{{dataset}}/eval/iGenVar_SL2/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
                               min_qual = min_qual_iGenVar),
        iGenVar_SL3   = expand("results/caller_comparison_short_read/{{dataset}}/eval/iGenVar_SL3/no_DUP_and_INV.min_qual_{min_qual}/pr_rec.txt",
                               min_qual = min_qual_iGenVar),
        Vaquita       = expand("results/caller_comparison_short_read/{{dataset}}/eval/Vaquita/min_qual_{min_qual}/pr_rec.txt",
                               min_qual=list(range(config["quality_ranges"]["from"],
                                                   config["quality_ranges"]["to"],
                                                   config["quality_ranges"]["step"]))),
        VaquitaLR_S   = expand("results/caller_comparison_short_read/{{dataset}}/eval/VaquitaLR_S/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                               min_qual = min_qual_VaquitaLR),
        VaquitaLR_SL1 = expand("results/caller_comparison_short_read/{{dataset}}/eval/VaquitaLR_SL1/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                               min_qual = min_qual_VaquitaLR),
        VaquitaLR_SL2 = expand("results/caller_comparison_short_read/{{dataset}}/eval/VaquitaLR_SL2/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                               min_qual = min_qual_VaquitaLR),
        VaquitaLR_SL3 = expand("results/caller_comparison_short_read/{{dataset}}/eval/VaquitaLR_SL3/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                               min_qual = min_qual_VaquitaLR),
        Delly2        = expand("results/caller_comparison_short_read/{{dataset}}/eval/Delly2/min_qual_{min_qual}/pr_rec.txt",
                               min_qual=list(range(config["quality_ranges"]["Delly2"]["from"],
                                                   config["quality_ranges"]["Delly2"]["to"],
                                                   config["quality_ranges"]["Delly2"]["step"]))),
        GRIDSS        = expand("results/caller_comparison_short_read/{{dataset}}/eval/GRIDSS/min_qual_{min_qual}/pr_rec.txt",
                               min_qual=list(range(config["quality_ranges"]["GRIDSS"]["from"],
                                                   config["quality_ranges"]["GRIDSS"]["to"],
                                                   config["quality_ranges"]["GRIDSS"]["step"]))),
        # Floating point exception
        # Generating GC wig file
        # Constructed GC wig in 189.84563970565796 sec
        # Loading GC wig file
        # Loading coverage wig file
        # Traceback (most recent call last):
        # File "/group/ag_abi/lbuntrock/anaconda3/envs/benchmarks/bin/tiddit", line 91, in <module>
        #     TIDDIT_calling.cluster(args)
        # File "TIDDIT_calling.py", line 207, in TIDDIT_calling.cluster
        # File "TIDDIT_coverage.py", line 8, in TIDDIT_coverage.coverage
        # FileNotFoundError: [Errno 2] No such file or directory: 'results/caller_comparison_short_read/hg38_Sim_default_2/TIDDIT/variants.wig'
        # Command exited with non-zero status 1
        TIDDIT        = expand("results/caller_comparison_short_read/{{dataset}}/eval/TIDDIT/min_qual_{min_qual}/pr_rec.txt",
                               min_qual=list(range(config["quality_ranges"]["TIDDIT"]["from"],
                                                   config["quality_ranges"]["TIDDIT"]["to"],
                                                   config["quality_ranges"]["TIDDIT"]["step"])))
    output:
        iGenVar_S     = temp("results/caller_comparison_short_read/{{dataset}}/eval/igenvar_s.all_results.txt"),
        iGenVar_SL1   = temp("results/caller_comparison_short_read/{{dataset}}/eval/igenvar_sl1.all_results.txt"),
        iGenVar_SL2   = temp("results/caller_comparison_short_read/{{dataset}}/eval/igenvar_sl2.all_results.txt"),
        iGenVar_SL3   = temp("results/caller_comparison_short_read/{{dataset}}/eval/igenvar_sl3.all_results.txt"),
        Vaquita       = temp("results/caller_comparison_short_read/{{dataset}}/eval/Vaquita.all_results.txt"),
        VaquitaLR_S   = temp("results/caller_comparison_short_read/{{dataset}}/eval/vaquitaLR_S.all_results.txt"),
        VaquitaLR_SL1 = temp("results/caller_comparison_short_read/{{dataset}}/eval/vaquitaLR_SL1.all_results.txt"),
        VaquitaLR_SL2 = temp("results/caller_comparison_short_read/{{dataset}}/eval/vaquitaLR_SL2.all_results.txt"),
        VaquitaLR_SL3 = temp("results/caller_comparison_short_read/{{dataset}}/eval/vaquitaLR_SL3.all_results.txt"),
        Delly2        = temp("results/caller_comparison_short_read/{{dataset}}/eval/Delly2.all_results.txt"),
        GRIDSS        = temp("results/caller_comparison_short_read/{{dataset}}/eval/GRIDSS.all_results.txt"),
        TIDDIT        = temp("results/caller_comparison_short_read/{{dataset}}/eval/TIDDIT.all_results.txt"),
        all = "results/caller_comparison_short_read/{dataset}/eval/all_results.txt"
    threads: 1
    run:
        shell("cat {input.iGenVar_S} > {output.iGenVar_S}")
        shell("cat {input.iGenVar_SL1} > {output.iGenVar_SL1}")
        shell("cat {input.iGenVar_SL2} > {output.iGenVar_SL2}")
        shell("cat {input.iGenVar_SL3} > {output.iGenVar_SL3}")
        shell("cat {input.Vaquita} > {output.Vaquita}")
        shell("cat {input.VaquitaLR_S} > {output.VaquitaLR_S}")
        shell("cat {input.VaquitaLR_SL1} > {output.VaquitaLR_SL1}")
        shell("cat {input.VaquitaLR_SL2} > {output.VaquitaLR_SL2}")
        shell("cat {input.VaquitaLR_SL3} > {output.VaquitaLR_SL3}")
        shell("cat {input.Delly2} > {output.Delly2}")
        shell("cat {input.GRIDSS} > {output.GRIDSS}")
        shell("cat {input.TIDDIT} > {output.TIDDIT}")
        shell("""
            cat {output.iGenVar_S} {output.iGenVar_SL1} {output.iGenVar_SL2} {output.iGenVar_SL3} \
                {output.Vaquita} {output.VaquitaLR_S} {output.VaquitaLR_SL1} {output.VaquitaLR_SL2} {output.VaquitaLR_SL3} \
                {output.Delly2} {output.GRIDSS} {output.TIDDIT} > {output.all}
        """)
