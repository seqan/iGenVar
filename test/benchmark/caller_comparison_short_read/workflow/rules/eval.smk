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
    log:
        "logs/caller_comparison_short_read/truvari/truvari_output.{dataset}.{caller}.{min_qual}.log"
    run:
        if wildcards.dataset == 'Illumina_Paired_End':
            truth_set_gz = config["truth_set"]["gz"],
            truth_set_bed = config["truth_set"]["bed"]
        else: # wildcards.dataset == 'Illumina_Mate_Pair'
            truth_set_gz = config["truth_set_renamed_chr"]["gz"],
            truth_set_bed = config["truth_set_renamed_chr"]["bed"]
        shell("""
            rm -rf {params.output_dir} && truvari bench -b {truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
                --passonly --includebed {truth_set_bed} &>> {log}
        """)
        # -f data/reference/hs37d5.fa

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
        # [W::vcf_parse_info] INFO 'CE' is not defined in the header, assuming Type=String
        # [W::bcf_update_info] INFO/END=0 is smaller than POS at chr1:1
        # ...
        # KeyError: 'unknown INFO: CE'
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
                               min_qual = min_qual_VaquitaLR)
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
        shell("""
            cat {output.iGenVar_S} {output.iGenVar_SL1} {output.iGenVar_SL2} {output.iGenVar_SL3} \
                {output.Vaquita} {output.VaquitaLR_S} {output.VaquitaLR_SL1} {output.VaquitaLR_SL2} {output.VaquitaLR_SL3} \
                > {output.all}
        """)
