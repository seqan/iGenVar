min_qual_iGenVar = list(range(config["quality_ranges"]["iGenVar"]["from"],
                              config["quality_ranges"]["iGenVar"]["to"],
                              config["quality_ranges"]["iGenVar"]["step"]))

rule DUP_as_INS:
    input:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.vcf"
    output:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.DUP_as_INS.vcf"
    shell:
        "sed -e 's/<DUP:TANDEM>/<INS>/g' {input.vcf} | sed -e 's/SVTYPE=DUP/SVTYPE=INS/g' > {output.vcf}"

rule DUP_as_INS_filter_vcf:
    input:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.DUP_as_INS.vcf"
    output:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.DUP_as_INS.min_qual_{min_qual}.vcf"
    conda:
        "../../../envs/bcftools.yaml"
    shell:
        "bcftools view -i 'QUAL>={wildcards.min_qual}.00' {input.vcf} > {output.vcf}"

rule DUP_as_INS_truvari:
    input:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.DUP_as_INS.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.DUP_as_INS.min_qual_{min_qual}.vcf.gz.tbi"
    output:
        summary = "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/DUP_as_INS.min_qual_{min_qual}/summary.txt"
    log:
        "logs/truvari/truvari_output.{dataset}.{long_read_enhancement}.{min_qual}.log"
    params:
        output_dir = "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/DUP_as_INS.min_qual_{min_qual}",
        truth_set_gz = config["truth_set_HG002_renamed_chr"]["gz"],
        truth_set_bed = config["truth_set_HG002_renamed_chr"]["bed"]
    shell:
        """
        rm -rf {params.output_dir} && truvari bench -b {params.truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
            --passonly --includebed {params.truth_set_bed} &>> {log}
        """
        # -f data/reference/hs37d5.fa

rule DUP_as_INS_reformat_truvari_results:
    input:
        "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/DUP_as_INS.min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.long_read_enhancement}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule DUP_as_INS_cat_truvari_results_all:
    input:
        S      = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/no_enhancement/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x1  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x1/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x2  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x2/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x3  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x3/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x5  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x5/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x10 = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x10/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x1   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x1/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x2   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x2/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x3   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x3/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x5   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x5/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x10  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x10/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2    = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/full_enhancement/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2     = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2/DUP_as_INS.min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar)
    output:
        S      = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.S.all_results.txt"),
        SL2x1  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.SL2x1.all_results.txt"),
        SL2x2  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.SL2x2.all_results.txt"),
        SL2x3  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.SL2x3.all_results.txt"),
        SL2x5  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.SL2x5.all_results.txt"),
        SL2x10 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.SL2x10.all_results.txt"),
        L2x1   = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.L2x1.all_results.txt"),
        L2x2   = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.L2x2.all_results.txt"),
        L2x3   = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.L2x3.all_results.txt"),
        L2x5   = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.L2x5.all_results.txt"),
        L2x10  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.L2x10.all_results.txt"),
        SL2    = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.SL2.all_results.txt"),
        L2     = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/DUP_as_INS.L2.all_results.txt"),
        all = "results/caller_comparison_long_read_enhancement/{dataset}/eval/DUP_as_INS.all_results.txt"
    run:
        shell("cat {input.S} > {output.S}")
        shell("cat {input.SL2x1} > {output.SL2x1}")
        shell("cat {input.SL2x2} > {output.SL2x2}")
        shell("cat {input.SL2x3} > {output.SL2x3}")
        shell("cat {input.SL2x5} > {output.SL2x5}")
        shell("cat {input.SL2x10} > {output.SL2x10}")
        shell("cat {input.L2x1} > {output.L2x1}")
        shell("cat {input.L2x2} > {output.L2x2}")
        shell("cat {input.L2x3} > {output.L2x3}")
        shell("cat {input.L2x5} > {output.L2x5}")
        shell("cat {input.L2x10} > {output.L2x10}")
        shell("cat {input.SL2} > {output.SL2}")
        shell("cat {input.L2} > {output.L2}")
        shell("""
            cat {output.S} {output.SL2x1} {output.SL2x2} {output.SL2x3} {output.SL2x5} {output.SL2x10} \
                {output.L2x1} {output.L2x2} {output.L2x3} {output.L2x5} {output.L2x10} \
                {output.SL2} {output.L2} > {output.all}
        """)
