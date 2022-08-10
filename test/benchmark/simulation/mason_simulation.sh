#!/usr/bin/env sh

cd data
mkdir -p simulation_GRCh38 && cd simulation_GRCh38

conda env create -f Repos/iGenVar/test/benchmark/envs/simulation.yaml
conda activate simulation

# With a 47644184 bp of chr21 and 10000000 simulated reads, which are 200 bp long, we have a coverage of:
# (10000000*200)/47644184 = 41.9778414087

# Cut chromosome 21 from reference
sed -n '/^>chr21/,$p' ../reference/GRCh38/hg38.fa > hg38_chr21_.fa
sed '/>chr22/,$d' hg38_chr21_.fa > hg38_chr21.fa

bwa index hg38_chr21.fa

# Mason
# Tools for Biological Sequence Simulation
# ========================================
# https://github.com/seqan/seqan/blob/master/apps/mason2/README

git clone --single-branch --depth 1 --branch develop https://github.com/seqan/seqan
cd ../.. && mkdir seqan && cd seqan
cmake -DCMAKE_BUILD_TYPE=Release ../../Repos/seqan
make -j16 -k mason_simulator mason_variator

# mason_variator - Variation Simulation
# =====================================
# mason_variator [OPTIONS] -ir IN.fa -ov OUT.vcf [-of OUT.fa]
#   Variation Simulation:
#     --snp-rate DOUBLE
#           Per-base SNP rate. In range [0.0..1.0]. Default: 0.0001.
#     --small-indel-rate DOUBLE
#           Small indel rate. In range [0.0..1.0]. Default: 0.000001.
#     --min-small-indel-size INTEGER
#           Minimal small indel size to simulate. In range [0..inf]. Default: 1.
#     --max-small-indel-size INTEGER
#           Maximal small indel size to simulate. In range [0..inf]. Default: 6.
#     --sv-indel-rate DOUBLE
#           Per-base SNP rate. In range [0.0..1.0]. Default: 0.0000001.
#     --sv-inversion-rate DOUBLE
#           Per-base SNP rate. In range [0.0..1.0]. Default: 0.0000001.
#     --sv-translocation-rate DOUBLE
#           Per-base SNP rate. In range [0.0..1.0]. Default: 0.0000001.
#     --sv-duplication-rate DOUBLE
#           Per-base SNP rate. In range [0.0..1.0]. Default: 0.0000001.
#     --min-sv-size INTEGER
#           Minimal SV size to simulate. In range [0..inf]. Default: 50.
#     --max-sv-size INTEGER
#           Maximal SV size to simulate. In range [0..inf]. Default: 1000.
SMALL_INDEL_SIZE=19
MIN_SV_SIZE=20
MAX_SV_SIZE=20000
SV_RATE=0.00001

./../../build/seqan/bin/mason_variator --in-reference hg38_chr21.fa --out-vcf hg38_SV_simulation_default.vcf --seed 0 \
                                       > hg38_SV_simulation_default.log
./../../build/seqan/bin/mason_variator --in-reference hg38_chr21.fa --out-vcf hg38_SV_simulation_InDel.vcf --seed 1 \
                                       --snp-rate 0 --small-indel-rate ${SV_RATE} --max-small-indel-size ${SMALL_INDEL_SIZE} \
                                       --sv-indel-rate ${SV_RATE} --min-sv-size ${MIN_SV_SIZE} --max-sv-size ${MAX_SV_SIZE} \
                                       --sv-inversion-rate 0 --sv-translocation-rate 0 --sv-duplication-rate 0 \
                                       > hg38_SV_simulation_InDel.log
./../../build/seqan/bin/mason_variator --in-reference hg38_chr21.fa --out-vcf hg38_SV_simulation_noSNP.vcf --seed 2 \
                                       --snp-rate 0 --small-indel-rate ${SV_RATE} --max-small-indel-size ${SMALL_INDEL_SIZE} \
                                       --sv-indel-rate ${SV_RATE} --min-sv-size ${MIN_SV_SIZE} --max-sv-size ${MAX_SV_SIZE} \
                                       --sv-inversion-rate ${SV_RATE} --sv-translocation-rate ${SV_RATE} --sv-duplication-rate ${SV_RATE} \
                                       > hg38_SV_simulation_noSNP.log
./../../build/seqan/bin/mason_variator --in-reference hg38_chr21.fa --out-vcf hg38_SV_simulation_SNPandSV.vcf --seed 3 \
                                       --snp-rate ${SV_RATE} --small-indel-rate ${SV_RATE} --max-small-indel-size ${SMALL_INDEL_SIZE} \
                                       --sv-indel-rate ${SV_RATE} --min-sv-size ${MIN_SV_SIZE} --max-sv-size ${MAX_SV_SIZE} \
                                       --sv-inversion-rate ${SV_RATE} --sv-translocation-rate ${SV_RATE} --sv-duplication-rate ${SV_RATE} \
                                       > hg38_SV_simulation_SNPandSV.log

# sort the vcf truth
picard SortVcf -I hg38_SV_simulation_default.vcf \
               -O hg38_SV_simulation_default_sorted.vcf -Xms1g -Xmx100g --TMP_DIR tmp/picard/ \
               >> hg38_SV_simulation_default.log
picard SortVcf -I hg38_SV_simulation_InDel.vcf \
               -O hg38_SV_simulation_InDel_sorted.vcf -Xms1g -Xmx100g --TMP_DIR tmp/picard/ \
               >> hg38_SV_simulation_InDel.log
picard SortVcf -I hg38_SV_simulation_noSNP.vcf \
               -O hg38_SV_simulation_noSNP_sorted.vcf -Xms1g -Xmx100g --TMP_DIR tmp/picard/ \
               >> hg38_SV_simulation_noSNP.log
picard SortVcf -I hg38_SV_simulation_SNPandSV.vcf \
               -O hg38_SV_simulation_SNPandSV_sorted.vcf -Xms1g -Xmx100g --TMP_DIR tmp/picard/ \
               >> hg38_SV_simulation_InDel.log

bgzip -c hg38_SV_simulation_default_sorted.vcf > hg38_SV_simulation_default_sorted.vcf.gz
bgzip -c hg38_SV_simulation_InDel_sorted.vcf > hg38_SV_simulation_InDel_sorted.vcf.gz
bgzip -c hg38_SV_simulation_noSNP_sorted.vcf > hg38_SV_simulation_noSNP_sorted.vcf.gz
bgzip -c hg38_SV_simulation_SNPandSV_sorted.vcf > hg38_SV_simulation_SNPandSV_sorted.vcf.gz

# create index
tabix -p vcf hg38_SV_simulation_default_sorted.vcf.gz
tabix -p vcf hg38_SV_simulation_InDel_sorted.vcf.gz
tabix -p vcf hg38_SV_simulation_noSNP_sorted.vcf.gz
tabix -p vcf hg38_SV_simulation_SNPandSV_sorted.vcf.gz

# mason_simulator - Read Simulation
# =================================
# mason_simulator [OPTIONS] -ir IN.fa --num-fragments NUM [-iv IN.vcf] --out LEFT.fq [--out-right RIGHT.fq]
# Simulate NUM reads/pairs from the reference sequence IN.fa, potentially with variants from IN.vcf. In case that
#       both -o and -or are given, write out paired-end data, if only -io is given, only single-end reads are simulated.
# -oa, --out-alignment OUTPUT_FILE
#       SAM/BAM file with alignments. Valid filetypes are: .sam[.*] and .bam, where * is any of the following
#       extensions: gz and bgzf for transparent (de)compression.
# -ir, --input-reference INPUT_FILE
#       Path to FASTA file to read the reference from. Valid filetypes are: .sam[.*], .raw[.*], .gbk[.*], .frn[.*],
#       .fq[.*], .fna[.*], .ffn[.*], .fastq[.*], .fasta[.*], .faa[.*], .fa[.*], .embl[.*], and .bam, where * is any
#       of the following extensions: gz and bgzf for transparent (de)compression.
# -iv, --input-vcf INPUT_FILE
#       Path to the VCF file with variants to apply. Valid filetype is: .vcf[.*], where * is any of the following
#       extensions: gz and bgzf for transparent (de)compression.

# ! needs indexed fasta file !
./../../build/seqan/bin/mason_simulator --input-reference hg38_chr21.fa --input-vcf hg38_SV_simulation_default.vcf \
                                        --out hg38_simulation_10000000_Illumina_reads_default.fastq \
                                        --num-fragments 10000000 --seq-technology illumina --illumina-read-length 200 \
                                        >> hg38_SV_simulation_default.log
./../../build/seqan/bin/mason_simulator --input-reference hg38_chr21.fa --input-vcf hg38_SV_simulation_InDel.vcf \
                                        --out hg38_simulation_10000000_Illumina_reads_InDel.fastq \
                                        --num-fragments 10000000 --seq-technology illumina --illumina-read-length 200 \
                                        >> hg38_SV_simulation_InDel.log
./../../build/seqan/bin/mason_simulator --input-reference hg38_chr21.fa --input-vcf hg38_SV_simulation_noSNP.vcf \
                                        --out hg38_simulation_10000000_Illumina_reads_noSNP.fastq \
                                        --num-fragments 10000000 --seq-technology illumina --illumina-read-length 200 \
                                        >> hg38_SV_simulation_noSNP.log
./../../build/seqan/bin/mason_simulator --input-reference hg38_chr21.fa --input-vcf hg38_SV_simulation_SNPandSV.vcf \
                                        --out hg38_simulation_10000000_Illumina_reads_SNPandSV.fastq \
                                        --num-fragments 10000000 --seq-technology illumina --illumina-read-length 200 \
                                        >> hg38_SV_simulation_InDel.log

bwa mem hg38_chr21.fa hg38_simulation_10000000_Illumina_reads_default.fastq \
                    > hg38_simulation_10000000_Illumina_reads_default.sam
bwa mem hg38_chr21.fa hg38_simulation_10000000_Illumina_reads_InDel.fastq \
                    > hg38_simulation_10000000_Illumina_reads_InDel.sam
bwa mem hg38_chr21.fa hg38_simulation_10000000_Illumina_reads_noSNP.fastq \
                    > hg38_simulation_10000000_Illumina_reads_noSNP.sam
bwa mem hg38_chr21.fa hg38_simulation_10000000_Illumina_reads_SNPandSV.fastq \
                    > hg38_simulation_10000000_Illumina_reads_SNPandSV.sam

# You can use the --mate-orientation to set the relative orientation when doing paired-end sequencing. The valid
# values are given in the following.

# FR    Reads are inward-facing, the same as Illumina paired-end reads: R1 --> <-- R2.
# RF    Reads are outward-facing, the same as Illumina mate-pair reads: R1 <-- --> R2.
# FF    Reads are on the same strand: R1 --> --> R2.
# FF2   Reads are on the same strand but the "right" reads are sequenced to the left of the "left" reads, same as
#         454 paired: R2 --> --> R1.

samtools sort -O BAM hg38_simulation_10000000_Illumina_reads_default.sam \
                   > hg38_simulation_10000000_Illumina_reads_default_sorted.bam
samtools sort -O BAM hg38_simulation_10000000_Illumina_reads_InDel.sam \
                   > hg38_simulation_10000000_Illumina_reads_InDel_sorted.bam
samtools sort -O BAM hg38_simulation_10000000_Illumina_reads_noSNP.sam \
                   > hg38_simulation_10000000_Illumina_reads_noSNP_sorted.bam
samtools sort -O BAM hg38_simulation_10000000_Illumina_reads_SNPandSV.sam \
                   > hg38_simulation_10000000_Illumina_reads_SNPandSV_sorted.bam

samtools index hg38_simulation_10000000_Illumina_reads_default_sorted.bam
samtools index hg38_simulation_10000000_Illumina_reads_InDel_sorted.bam
samtools index hg38_simulation_10000000_Illumina_reads_noSNP_sorted.bam
samtools index hg38_simulation_10000000_Illumina_reads_SNPandSV_sorted.bam

cd ../..
