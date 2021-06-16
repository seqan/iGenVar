localrules: bgzip, tabix, callset_eval_svim, callset_eval, reformat_truvari_results, reformat_truvari_results_svim, cat_truvari_results_all, cat_truvari_results_full, cat_truvari_results_svim_parameters

def get_vcf(wildcards):
    return config["truth"][wildcards.vcf.split(".")[0]]

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


rule callset_eval_svim:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/SVIM/{parameters}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/SVIM/{parameters}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        summary="pipeline/SVIM_results/{parameters}/{minscore}/{vcf}/summary.txt"
    params:
        out_dir="pipeline/SVIM_results/{parameters}/{minscore}/{vcf}",
        vcf="{vcf}"
    threads: 1
    log:
        log="logs/truvari/truvari.svim.{parameters}.{minscore}.{vcf}.log"
    conda:
        "../../../envs/truvari.yaml"
    script:
        "../scripts/run_truvari.py"


rule callset_eval:
    input:
        genome = config["reference"],
        truth_vcf = get_vcf,
        truth_bed = config["truth"]["bed"],
        calls = "pipeline/{caller}/min_{minscore}.indel.vcf.gz",
        index = "pipeline/{caller}/min_{minscore}.indel.vcf.gz.tbi"
    output:
        summary="pipeline/{caller,Sniffles|pbsv}_results/{minscore}/{vcf}/summary.txt"
    params:
        out_dir="pipeline/{caller}_results/{minscore}/{vcf}",
        vcf="{vcf}"
    threads: 1
    log:
        log="logs/truvari/truvari.{caller}.{minscore}.{vcf}.log"
    conda:
        "../../../envs/truvari.yaml"
    script:
        "../scripts/run_truvari.py"


rule reformat_truvari_results:
    input:
        "pipeline/{caller}_results/{minscore}/{vcf}/summary.txt"
    output:
        "pipeline/{caller,Sniffles|pbsv}_results/{minscore}/{vcf}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"{wildcards.caller}\", \"{wildcards.data}\", \"{wildcards.vcf}\", {wildcards.minscore}, $1, $2 }}' > {output}"

rule reformat_truvari_results_svim:
    input:
        "pipeline/SVIM_results/{parameters}/{minscore}/{vcf}/summary.txt"
    output:
        "pipeline/SVIM_results/{parameters}/{minscore}/{vcf}/pr_rec.txt"
    threads: 1
    shell:
        "cat {input} | grep 'precision\|recall' | tr -d ',' |sed 's/^[ \t]*//' | tr -d '\"' | tr -d ' ' | tr ':' '\t' | awk 'OFS=\"\\t\" {{ print \"SVIM\", \"{wildcards.data}\", \"{wildcards.vcf}\", {wildcards.minscore}, $1, $2 }}' > {output}"


rule cat_truvari_results_all:
    input:
        svim = expand("pipeline/SVIM_results/1000_900_1.0_0.5/{minscore}/{vcf}/pr_rec.txt",
                          minscore=[0] + SVIM_THRESHOLDS, vcf=VCFS),
        sniffles = expand("pipeline/Sniffles_results/{minscore}/{vcf}/pr_rec.txt",
                          minscore=list(range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"]+1, config["minimums"]["sniffles_step"])),
                          vcf=VCFS),
        pbsv = expand("pipeline/pbsv_results/{minscore}/{vcf}/pr_rec.txt",
                          minscore=list(range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"]+1, config["minimums"]["pbsv_step"])),
                          vcf=VCFS)
    output:
        svim = temp("pipeline/eval/svim.all_results.txt"),
        sniffles = temp("pipeline/eval/sniffles.all_results.txt"),
        pbsv = temp("pipeline/eval/pbsv.all_results.txt"),
        all = "pipeline/eval/all_results.txt"
    threads: 1
    run:
        shell("cat {input.svim} > {output.svim}")
        shell("cat {input.sniffles} > {output.sniffles}")
        shell("cat {input.pbsv} > {output.pbsv}")
        shell("cat {output.svim} {output.sniffles} {output.pbsv} > {output.all}")

rule cat_truvari_results_full:
    input:
        svim = expand("pipeline/SVIM_results/pooled/1000_900_1.0_0.5/{minscore}/{vcf}/pr_rec.txt",
                          minscore=[0] + SVIM_THRESHOLDS, vcf=VCFS),
        sniffles = expand("pipeline/Sniffles_results/pooled/{minscore}/{vcf}/pr_rec.txt",
                          minscore=list(range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"]+1, config["minimums"]["sniffles_step"])),
                          vcf=VCFS),
        pbsv = expand("pipeline/pbsv_results/pooled/{minscore}/{vcf}/pr_rec.txt",
                          minscore=list(range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"]+1, config["minimums"]["pbsv_step"])),
                          vcf=VCFS)
    output:
        svim = temp("pipeline/eval/svim.full_results.txt"),
        sniffles = temp("pipeline/eval/sniffles.full_results.txt"),
        pbsv = temp("pipeline/eval/pbsv.full_results.txt"),
        all = "pipeline/eval/full_results.txt"
    threads: 1
    run:
        shell("cat {input.svim} > {output.svim}")
        shell("cat {input.sniffles} > {output.sniffles}")
        shell("cat {input.pbsv} > {output.pbsv}")
        shell("cat {output.svim} {output.sniffles} {output.pbsv} > {output.all}")

rule cat_truvari_results_svim_parameters:
    input:
        svim = expand("pipeline/SVIM_results/{pmd}_{pdn}_{edn}_{cmd}/{minscore}/{vcf}/pr_rec.txt",
                          pmd = [1000],
                          pdn = [900],
                          edn = [1.0, 1.5, 2.0, 2.5, 3.0, 4.0],
                          cmd = [0.2, 0.25, 0.3, 0.35, 0.4, 0.45, 0.5],
                          minscore= SVIM_THRESHOLDS,
                          vcf=VCFS)
    output:
        all = "pipeline/eval/svim_parameter_results.txt"
    threads: 1
    run:
        with open(output.all, 'w') as output_file:
            for f in input.svim:
                parameters = f.split("/")[4]
                with open(f, 'r') as input_file:
                    for line in input_file:
                        print("%s\t%s" % (parameters, line.strip()), file=output_file)
