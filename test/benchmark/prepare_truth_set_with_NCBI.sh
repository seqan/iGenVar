TRUTH_SET=HG002_SVs_Tier1_v0.6

# NCBI Genome Remapping Service

conda env create -f Repos/iGenVar/test/benchmark/envs/perl.yaml

conda create -n perl -c conda-forge -c bioconda perl perl-xml-xpath perl-test-leaktrace perl-html-parser perl-io-socket-ssl # This is perl 5, version 32, subversion 1 (v5.32.1)
conda activate perl

# download mapping script (perl)
${WGET} ${NCBI}/pub/remap/remap_api.pl
# get required perl modules:
# Getopt::Long
cpan -i Getopt::Long           # Getopt::Long is up to date (2.52).
# LWP::UserAgent
cpan -i LWP::UserAgent      # LWP::UserAgent is up to date (6.62).
cpan -i LWP::Protocol::https
# HTTP::Request::Common qw(POST), qw(GET)
cpan -i HTTP::Request::Common  # HTTP::Request::Common is up to date (6.36).
# HTTP::Headers
cpan -i HTTP::Headers          # HTTP::Headers is up to date (6.36).
# XML::XPath
cpan -i YAML
cpan -i XML::Parser
cpan -i XML::XPath             # XML::XPath does not compile
# XML::XPath::XMLPars
cpan -i XML::XPath::XMLPars # Skipping XML::XPath::XMLPars because I couldn't find a matching namespace.

## split vcf
conda install -c bioconda vcftools
vcftools --vcf truth_set/${TRUTH_SET}.vcf --chr 1 --chr 2 --chr 3 --chr 4 --chr 5 --recode --recode-INFO-all --out Chr1-2-3-4-5
vcftools --vcf truth_set/${TRUTH_SET}.vcf --chr 6 --chr 7 --chr 8 --chr 9 --chr 10 --chr 11 --chr 12 --recode --recode-INFO-all --out Chr6-7-8-9-10-11-12
vcftools --vcf truth_set/${TRUTH_SET}.vcf --chr 13 --chr 14 --chr 15 --chr 16 --chr 17 --chr 18 --chr 19 --chr 20 --chr 21 --chr 22 --chr X --chr Y --chr MT --recode --recode-INFO-all --out Chr13-14-15-16-17-18-19-20-21-22-X-Y-MT

bcftools merge


# source
${WGET} ${NCBI}/pub/remap/Homo_sapiens/2.2/GCF_000001405.13_GRCh37/GCF_000001405.39_GRCh38.p13/GCF_000001405.13-GCF_000001405.39.gff
# target
${WGET} ${NCBI}/pub/remap/Homo_sapiens/2.2/GCF_000001405.26_GRCh38/GCF_000001405.39_GRCh38.p13/GCF_000001405.26-GCF_000001405.39.gff

echo "$(tput setaf 1)$(tput setab 7)------- NCBI Genome Remapping Service installed and prepared (5.4/5.6) --------$(tput sgr 0)" 1>&3

# convert truth set (VCF & BED)

# --mode asm-asm: remap mode Assembly-Assembly
perl remap_api.pl --mode asm-asm --from GCF_000001405.13 --dest GCF_000001405.39 --annotation ${TRUTH_SET}.vcf --annot_out ${TRUTH_SET}.Hg38.vcf --report_out ${TRUTH_SET}_GRCh37_GRCh38.tsv --gbench_out ${TRUTH_SET}_GRCh38.gbp
perl remap_api.pl --mode asm-asm --from GCF_000001405.13 --dest GCF_000001405.39 --annotation ${TRUTH_SET}.bed --annot_out ${TRUTH_SET}.Hg38.bed --report_out ${TRUTH_SET}_GRCh37_GRCh38.tsv --gbench_out ${TRUTH_SET}_GRCh38.gbp

echo "$(tput setaf 1)$(tput setab 7)------- truth set files remapped (5.5/5.6) --------$(tput sgr 0)" 1>&3

## VCF
# create new tbi file
# tabix -p vcf truth_set/${TRUTH_SET}.renamed_chr.Hg38.vcf.gz
# unzip vcf
bgzip -d -c truth_set/${TRUTH_SET}.Hg38.vcf.gz > truth_set/${TRUTH_SET}.Hg38.vcf
# bgzip -d -c truth_set/NA24385.GRCh38.large_svs.vcf.gz > truth_set/NA24385.GRCh38.large_svs.vcf

# add chr to chromosome names
sed -e 's/##contig=<ID=/##contig=<ID=chr/g' truth_set/${TRUTH_SET}.vcf > truth_set/${TRUTH_SET}.renamed_contig.vcf
awk '{if($0q!~ /^#/) print "chr"$0; else print $0}' truth_set/${TRUTH_SET}.renamed_contig.vcf > truth_set/${TRUTH_SET}.renamed_chr.vcf
sed -e 's/##contig=<ID=/##contig=<ID=chr/g' truth_set/${TRUTH_SET}.Hg38.vcf > truth_set/${TRUTH_SET}.Hg38.renamed_contig.vcf
awk '{if($0q!~ /^#/) print "chr"$0; else print $0}' truth_set/${TRUTH_SET}.Hg38.renamed_contig.vcf > truth_set/${TRUTH_SET}.Hg38.renamed_chr.vcf
# zip file again
bgzip -c truth_set/${TRUTH_SET}.Hg38.renamed_chr.vcf > truth_set/${TRUTH_SET}.Hg38.renamed_chr.vcf.gz
# create tbi files
tabix -p vcf truth_set/${TRUTH_SET}.Hg38.vcf.gz
tabix -p vcf truth_set/${TRUTH_SET}.Hg38.renamed_chr.vcf.gz
# tabix -p vcf truth_set/NA24385.GRCh38.large_svs.vcf.gz

# convert to HG38
CrossMap.py bed hg19ToHg38.over.chain.gz truth_set/${TRUTH_SET}.renamed_chr.bed truth_set/${TRUTH_SET}.Hg38.renamed_chr.bed
# remove chr from chromosome names again
sed -e 's!chr!!' truth_set/${TRUTH_SET}.renamed_chr.Hg38.bed > truth_set/${TRUTH_SET}.Hg38.bed

# convert2bed --input=vcf --insertions < NA24385.GRCh38.large_svs.vcf.gz > insertions.bed

echo "$(tput setaf 1)$(tput setab 7)------- truth set files prepared (5.6/5.6) --------$(tput sgr 0)" 1>&3
