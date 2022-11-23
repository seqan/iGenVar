min_qual_iGenVar = list(range(config["quality_ranges"]["iGenVar"]["from"],
                              config["quality_ranges"]["iGenVar"]["to"],
                              config["quality_ranges"]["iGenVar"]["step"]))

rule filter_vcf:
    input:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.vcf"
    output:
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.min_qual_{min_qual}.vcf"
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
        vcf = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.min_qual_{min_qual}.vcf.gz",
        index = "results/caller_comparison_long_read_enhancement/{dataset}/{long_read_enhancement}/variants.min_qual_{min_qual}.vcf.gz.tbi"
    output:
        summary = "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/min_qual_{min_qual}/summary.txt"
    log:
        "logs/truvari/truvari_output.{dataset}.{long_read_enhancement}.{min_qual}.log"
    params:
        output_dir = "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/min_qual_{min_qual}",
        truth_set_gz = config["truth_set_HG002_renamed_chr"]["gz"],
        truth_set_bed = config["truth_set_HG002_renamed_chr"]["bed"]
    shell:
        """
        rm -rf {params.output_dir} && truvari bench -b {params.truth_set_gz} -c {input.vcf} -o {params.output_dir} -p 0 \
            --passonly --includebed {params.truth_set_bed} &>> {log}
        """
        # -f data/reference/hs37d5.fa

rule reformat_truvari_results:
    input:
        "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/min_qual_{min_qual}/summary.txt"
    output:
        "results/caller_comparison_long_read_enhancement/{dataset}/eval/{long_read_enhancement}/min_qual_{min_qual}/pr_rec.txt"
    threads: 1
    shell:
        """
        cat {input} | grep '\<precision\>\|\<recall\>' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' \
            | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.long_read_enhancement}\", \"{wildcards.min_qual}\", $1, $2 }}' \
            > {output}
        """

rule cat_truvari_results_all:
    input:
        S      = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/no_enhancement/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x1  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x1/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x2  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x2/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x3  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x3/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x5  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x5/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2x10 = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x10/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x1   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x1/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x2   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x2/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x3   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x3/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x5   = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x5/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2x10  = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x10/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        SL2    = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/full_enhancement/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar),
        L2     = expand("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2/min_qual_{min_qual}/pr_rec.txt",
                        min_qual = min_qual_iGenVar)
    output:
        S      = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/S.all_results.txt"),
        SL2x1  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x1.all_results.txt"),
        SL2x2  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x2.all_results.txt"),
        SL2x3  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x3.all_results.txt"),
        SL2x5  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x5.all_results.txt"),
        SL2x10 = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2x10.all_results.txt"),
        L2x1   = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x1.all_results.txt"),
        L2x2   = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x2.all_results.txt"),
        L2x3   = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x3.all_results.txt"),
        L2x5   = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x5.all_results.txt"),
        L2x10  = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2x10.all_results.txt"),
        SL2    = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/SL2.all_results.txt"),
        L2     = temp("results/caller_comparison_long_read_enhancement/{{dataset}}/eval/L2.all_results.txt"),
        all = "results/caller_comparison_long_read_enhancement/{dataset}/eval/all_results.txt"
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
                {output.SL2} {output.L2}> {output.all}
        """)
