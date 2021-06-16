localrules: filter_svim, fix_sniffles, filter_insertions_and_deletions, filter_insertions_and_deletions_svim

# SVIM
rule run_svim:
    input:
        bam = config["long_bam"],
        bai = config["long_bai"],
        genome = config["reference"]
    output:
        "pipeline/SVIM/{pmd,[0-9]+}_{pdn,[0-9]+}_{edn,[0-9\.]+}_{cmd,[0-9\.]+}/variants.vcf"
    resources:
        mem_mb = 20000,
        time_min = 600,
        io_gb = 100
    params:
        working_dir = "pipeline/SVIM/{pmd}_{pdn}_{edn}_{cmd}/",
        min_sv_size = config["parameters"]["min_sv_size"]
    threads: 1
    conda:
        "../../../envs/svim.yaml"
    shell:
        "svim alignment --sample {wildcards.data} \
         --partition_max_distance {wildcards.pmd} \
         --position_distance_normalizer {wildcards.pdn} \
         --edit_distance_normalizer {wildcards.edn} \
         --cluster_max_distance {wildcards.cmd} \
         --min_sv_size {params.min_sv_size} \
         --segment_gap_tolerance 20 \
         --segment_overlap_tolerance 20 \
         --interspersed_duplications_as_insertions \
         --tandem_duplications_as_insertions \
         --read_names \
         --max_sv_size 1000000 \
         --verbose \
         {params.working_dir} {input.bam} {input.genome}"

rule filter_svim:
    input:
        "pipeline/SVIM/{parameters}/variants.vcf"
    output:
        temp("pipeline/SVIM/{parameters}/min_{minscore,[0-9]+}.vcf")
    threads: 1
    shell:
        "bcftools view -e \"GT=='0/0'\" {input} | \
         awk 'OFS=\"\\t\" {{ if($1 ~ /^#/) {{ print $0 }} \
         else {{ if($6>={wildcards.minscore}) {{ print $1, $2, $3, $4, $5, $6, \"PASS\", $8, $9, $10 }} }} }}' > {output}"

# SNIFFLES
rule run_sniffles:
    input:
        bam = config["long_bam"],
        bai = config["long_bai"]
    output:
        expand("pipeline/Sniffles/raw_{minsupport}.vcf",
                minsupport=list(range(config["minimums"]["sniffles_from"], config["minimums"]["sniffles_to"]+1, config["minimums"]["sniffles_step"])))
    resources:
        mem_mb = 400000,
        time_min = 1200,
        io_gb = 100
    params:
        min_sv_size = config["parameters"]["min_sv_size"],
        tmpdir = "300",
        sniffles_from = config["minimums"]["sniffles_from"],
        sniffles_to = config["minimums"]["sniffles_to"],
        sniffles_step = config["minimums"]["sniffles_step"],
        outdir = "pipeline/Sniffles/"
    threads: 30
    conda:
        "../../../envs/sniffles.yaml"
    shell:
        "bash workflow/scripts/run_sniffles.sh {input.bam} {input.bai} {params.sniffles_from} {params.sniffles_to} {params.sniffles_step} {params.min_sv_size} {threads} {params.outdir}"

#see https://github.com/spiralgenetics/truvari/issues/43
rule fix_sniffles:
    input:
        "pipeline/Sniffles/raw_{support}.vcf"
    output:
        "pipeline/Sniffles/min_{support,[0-9]+}.vcf"
    shell:
        "sed 's/##INFO=<ID=SUPTYPE,Number=A/##INFO=<ID=SUPTYPE,Number=./' {input} > {output}"

#PBSV
rule run_pbsv:
    input:
        bam = config["long_bam"],
        bai = config["long_bai"]
        genome = config["reference"],
    output:
        expand("pipeline/pbsv/min_{minsupport}.vcf",
                minsupport=list(range(config["minimums"]["pbsv_from"], config["minimums"]["pbsv_to"]+1, config["minimums"]["pbsv_step"])))
    resources:
        mem_mb = 400000,
        time_min = 2000,
        io_gb = 100
    params:
        min_sv_size = config["parameters"]["min_sv_size"],
        tmpdir = "1000",
        pbsv_from = config["minimums"]["pbsv_from"],
        pbsv_to = config["minimums"]["pbsv_to"],
        pbsv_step = config["minimums"]["pbsv_step"],
        outdir = "pipeline/pbsv/"
    threads: 15
    conda:
        "../../../envs/pbsv.yaml"
    shell:
        "bash workflow/scripts/run_pbsv.sh {input.bam} {input.bai} {input.genome} {params.pbsv_from} {params.pbsv_to} {params.pbsv_step} {params.min_sv_size} {threads} {params.outdir}"

#Split to SV classes
rule filter_insertions_and_deletions:
    input:
        "pipeline/{caller}/min_{minscore}.vcf"
    output:
        "pipeline/{caller,Sniffles|pbsv}/min_{minscore,[0-9]+}.indel.vcf"
    threads: 1
    shell:
        "bcftools view -i 'SVTYPE=\"DEL\" | SVTYPE=\"INS\"' {input} | bcftools sort > {output}"

rule filter_insertions_and_deletions_svim:
    input:
        "pipeline/SVIM/{parameters}/min_{minscore}.vcf"
    output:
        "pipeline/SVIM/{parameters}/min_{minscore,[0-9]+}.indel.vcf"
    threads: 1
    shell:
        "bcftools view -i 'SVTYPE=\"DEL\" | SVTYPE=\"INS\"' {input} | bcftools sort > {output}"
