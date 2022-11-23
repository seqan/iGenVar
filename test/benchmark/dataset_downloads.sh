WGET="wget --retry-connrefused --waitretry=30 --read-timeout=30 --timeout=30 --tries=20 --no-clobber --no-verbose"
NCBI="ftp://ftp.ncbi.nlm.nih.gov/" # "https://ftp-trace.ncbi.nlm.nih.gov/" "https://ftp-trace.ncbi.nlm.nih.gov/giab/ftp/data/"
HG002="AshkenazimTrio/HG002_NA24385_son"
GENOMES="genomes/all/GCA/000/001/405" # README: https://ftp-trace.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/README.txt
ReferenceSamples="ReferenceSamples/giab/data"

# -------- -------- get data and unzip -------- -------- #
mkdir -p short_reads && cd short_reads

# HG002 - NA24385:
mkdir -p GRCh37 && cd GRCh37
mkdir -p HG002 && cd HG002
## NIST Illumina 2x250bps Paired-end
### reference: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
### README: .../NIST_Illumina_2x250bps/novoalign_bams/README_update_feb2019
### coverage: ~40-50x (samtools depth average = 65.8713)
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_Illumina_2x250bps/novoalign_bams/HG002.hs37d5.2x250.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_Illumina_2x250bps/novoalign_bams/HG002.hs37d5.2x250.bam.bai
cd ../.. && mkdir -p GRCh38 && cd GRCh38
mkdir -p HG002 && cd HG002
### reference: GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
### README: .../NIST_Illumina_2x250bps/novoalign_bams/README_update_feb2019
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_Illumina_2x250bps/novoalign_bams/HG002.GRCh38.2x250.bam.bai

## NIST Stanford Illumina 6kb matepair
### Reads were then mapped to the hg19 reference genome from ucsc or the GRCh38 reference genome with decoy but no alts
### using bwa mem (Li 2013) with default settings, and duplicates were marked using samblaster (Faust 2014).
cd ../.. && cd GRCh37/HG002
### reference: hg19 from ucsc
### README: .../NIST_Stanford_Illumina_6kb_matepair/README.NIST_Stanford_Illumina_6kb_matepair
### coverage: 13-17x (samtools depth average = 14.1621)
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_Stanford_Illumina_6kb_matepair/bams/HG002.mate_pair.sorted.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_Stanford_Illumina_6kb_matepair/bams/HG002.mate_pair.sorted.bam.bai
cd ../.. && cd GRCh38/HG002
### reference: GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna
### README: .../NIST_Stanford_Illumina_6kb_matepair/README.NIST_Stanford_Illumina_6kb_matepair
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_Stanford_Illumina_6kb_matepair/bams/GRCh38/HG002.sorted.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_Stanford_Illumina_6kb_matepair/bams/GRCh38/HG002.sorted.bam.bai

## NIST HiSeq - NHGRI Illumina 300X
cd ../.. && cd GRCh37/HG002
### README: .../NIST_HiSeq_HG002_Homogeneity-10953946/README_NIST_Illumina_pairedend_HG002.txt
### coverage: 60x (samtools depth average = 59.0516)
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.hs37d5.60x.1.bam.bai
cd ../.. && cd GRCh38/HG002
### README: .../NIST_HiSeq_HG002_Homogeneity-10953946/README_NIST_Illumina_pairedend_HG002.txt
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}NIST_HiSeq_HG002_Homogeneity-10953946/NHGRI_Illumina300X_AJtrio_novoalign_bams/HG002.GRCh38.60x.1.bam.bai

cd ..

mkdir -p NA12878 && cd NA12878

## NA12878 WGS GRCh38 - Platinum Genomes
### Project HiSeq 2000: TruSeq PCR-Free (Platinum Genomes)
### coverage: 50x (samtools depth average = 47.6849)
curl 'https://basespace-data-east.s3-external-1.amazonaws.com/47c60b82082847db8cf4a7d900e530b4/NA12878_S1.bam?AWSAccessKeyId=AKIARPYQJSWQRDJKAWUT&Expires=1665567801&response-content-disposition=attachment%3Bfilename%3DNA12878_S1.bam&response-content-type=application%2Foctet-stream&Signature=ca6d7EvQOvop49gEh%2Fh6pTs2HIQ%3D' -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:102.0) Gecko/20100101 Firefox/102.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' -H 'Accept-Language: de,en-US;q=0.7,en;q=0.3' -H 'Accept-Encoding: gzip, deflate, br' -H 'Referer: https://basespace.illumina.com/' -H 'Connection: keep-alive' -H 'Upgrade-Insecure-Requests: 1' -H 'Sec-Fetch-Dest: document' -H 'Sec-Fetch-Mode: navigate' -H 'Sec-Fetch-Site: cross-site' -H 'Sec-Fetch-User: ?1' \
    --output NA12878_S1.bam
curl 'https://basespace-data-east.s3-external-1.amazonaws.com/47c60b82082847db8cf4a7d900e530b4/NA12878_S1.bam.bai?AWSAccessKeyId=AKIARPYQJSWQRDJKAWUT&Expires=1665568255&response-content-disposition=attachment%3Bfilename%3DNA12878_S1.bam.bai&response-content-type=application%2Foctet-stream&Signature=PQNrr3brzcLKkqZRgh1zw3AJcuM%3D' -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:102.0) Gecko/20100101 Firefox/102.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' -H 'Accept-Language: de,en-US;q=0.7,en;q=0.3' -H 'Accept-Encoding: gzip, deflate, br' -H 'Referer: https://basespace.illumina.com/' -H 'Connection: keep-alive' -H 'Upgrade-Insecure-Requests: 1' -H 'Sec-Fetch-Dest: document' -H 'Sec-Fetch-Mode: navigate' -H 'Sec-Fetch-Site: cross-site' -H 'Sec-Fetch-User: ?1' \
--output NA12878_S1.bam.bai

cd ../../..

echo "$(tput setaf 1)$(tput setab 7)------- short reads downloaded (3.1/3.4) --------$(tput sgr 0)" 1>&3

mkdir -p long_reads && cd long_reads
# HG002 - NA24385:
## MtSinai PacBio (minimap2)
mkdir -p GRCh37 && cd GRCh37
### reference: ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
### README: .../PacBio_MtSinai_NIST/PacBio_minimap2_bam/README.txt
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/PacBio_MtSinai_NIST/PacBio_minimap2_bam/HG002_PacBio_GRCh37.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/PacBio_MtSinai_NIST/PacBio_minimap2_bam/HG002_PacBio_GRCh37.bam.bai
cd .. && mkdir -p GRCh38 && cd GRCh38
### reference: ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
### README: .../PacBio_MtSinai_NIST/PacBio_minimap2_bam/README.txt
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/PacBio_MtSinai_NIST/PacBio_minimap2_bam/HG002_PacBio_GRCh38.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/PacBio_MtSinai_NIST/PacBio_minimap2_bam/HG002_PacBio_GRCh38.bam.bai

## PacBio CCS 10kb (minimap2)
cd .. && cd GRCh37
### reference: The reads were aligned to the hs37d5 reference
### README: .../PacBio_CCS_10kb/alignment/README.txt
### coverage: 30x (samtools depth average = 28.6346)
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/PacBio_CCS_10kb/alignment/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/PacBio_CCS_10kb/alignment/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam.bai
cd .. && cd GRCh38
### reference: reads were aligned to the GRCh38_no_alt_analysis reference using pbmm2
### README: .../PacBio_CCS_10kb/GRCh38_no_alt_analysis/README.txt
### coverage: 30x (samtools depth average = 28.8953)
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/PacBio_CCS_10kb/GRCh38_no_alt_analysis/HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/PacBio_CCS_10kb/GRCh38_no_alt_analysis/HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam.bai

## 10X Genomics
cd .. && cd GRCh37
### reference:
### README: .../10XGenomics/README_10Xgenomes.txt
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/10XGenomics/NA24385_phased_possorted_bam.bam
${WGET} ${NCBI}/${ReferenceSamples}/${HG002}/10XGenomics/NA24385_phased_possorted_bam.bam.bai

# HG005 - NA24631:
## MtSinai PacBio (minimap2)
# ${WGET} ${NCBI}/${ReferenceSamples}/ChineseTrio/HG005_NA24631_son/MtSinai_PacBio/PacBio_minimap2_bam/HG005_PacBio_GRCh38.bam

# ${WGET} ${NCBI}/${ReferenceSamples}/ChineseTrio/HG005_NA24631_son/MtSinai_PacBio/PacBio_minimap2_bam/HG005_PacBio_GRCh38.bam.bai

# ## PacBio SequelII CCS 11kb
# ${WGET} ${NCBI}/${ReferenceSamples}/ChineseTrio/HG005_NA24631_son/PacBio_SequelII_CCS_11kb/HG005.SequelII.pbmm2.hs37d5.whatshap.haplotag.10x.bam

# ${WGET} ${NCBI}/${ReferenceSamples}/ChineseTrio/HG005_NA24631_son/PacBio_SequelII_CCS_11kb/HG005.SequelII.pbmm2.hs37d5.whatshap.haplotag.10x.bam.bai

cd ../..

echo "$(tput setaf 1)$(tput setab 7)------- long reads downloaded (3.2/3.4) --------$(tput sgr 0)" 1>&3

mkdir -p reference && cd reference
mkdir -p GRCh37 && cd GRCh37

# 1000genomes Reference Genome Sequence (hs37d5), based on NCBI GRCh37
${WGET} ${NCBI}/1000genomes/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
# Alternative location ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz
gzip --decompress --keep hs37d5.fa.gz

# ucsc
${WGET} http://hgdownload.cse.ucsc.edu/goldenpath/hg19/bigZips/hg19.fa.gz
gzip --decompress --keep hg19.fa.gz

cd .. && mkdir -p GRCh38 && cd GRCh38
# README: https://ftp-trace.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/README_analysis_sets.txt
${WGET} ${NCBI}/${GENOMES}/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
gzip --decompress --keep GCA_000001405.15_GRCh38_no_alt_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_analysis_set.fa
${WGET} ${NCBI}/${GENOMES}/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
gzip --decompress --keep GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
mv GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fa

# GRCh38_reference_genome
#     https://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/GRCh38_reference_genome/GRCh38_full_analysis_set_plus_decoy_hla.fa

cd ../..

echo "$(tput setaf 1)$(tput setab 7)------- reference downloaded (3.3/3.4) --------$(tput sgr 0)" 1>&3

# -------- -------- get truth set ../../data -------- -------- #
mkdir -p truth_set && cd truth_set
mkdir -p HG002 && cd HG002

# HG002 - NA24385:
## GRCh37
${WGET} ${NCBI}/${ReferenceSamples}/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz
${WGET} ${NCBI}/${ReferenceSamples}/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.vcf.gz.tbi
${WGET} ${NCBI}/${ReferenceSamples}/AshkenazimTrio/analysis/NIST_SVs_Integration_v0.6/HG002_SVs_Tier1_v0.6.bed

## GRCh38
# ${WGET} ${NCBI}/${ReferenceSamples}/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh38/NA24385_300G/NA24385.GRCh38.large_svs.vcf.gz
# ${WGET} ${NCBI}/${ReferenceSamples}/AshkenazimTrio/analysis/10XGenomics_ChromiumGenome_LongRanger2.2_Supernova2.0.1_04122018/GRCh38/NA24385_300G/NA24385.GRCh38.large_svs.vcf.gz.tbi

cd ..
mkdir -p NA12878 && cd NA12878

## NA12878 WGS GRCh38 - Platinum Genomes
### Project HiSeq 2000: TruSeq PCR-Free (Platinum Genomes)
curl 'https://basespace-data-east.s3-external-1.amazonaws.com/47c60b82082847db8cf4a7d900e530b4/NA12878_S1.vcf.gz?AWSAccessKeyId=AKIARPYQJSWQRDJKAWUT&Expires=1665568301&response-content-disposition=attachment%3Bfilename%3DNA12878_S1.vcf.gz&response-content-type=application%2Fx-gzip&Signature=yF%2Fk2li9bNaCdbQxqoX0BKtIdww%3D' -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:102.0) Gecko/20100101 Firefox/102.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' -H 'Accept-Language: de,en-US;q=0.7,en;q=0.3' -H 'Accept-Encoding: gzip, deflate, br' -H 'Referer: https://basespace.illumina.com/' -H 'Connection: keep-alive' -H 'Upgrade-Insecure-Requests: 1' -H 'Sec-Fetch-Dest: document' -H 'Sec-Fetch-Mode: navigate' -H 'Sec-Fetch-Site: cross-site' -H 'Sec-Fetch-User: ?1' \
--output NA12878_S1.vcf.gz
curl 'https://basespace-data-east.s3-external-1.amazonaws.com/47c60b82082847db8cf4a7d900e530b4/NA12878_S1.vcf.gz.tbi?AWSAccessKeyId=AKIARPYQJSWQRDJKAWUT&Expires=1665568395&response-content-disposition=attachment%3Bfilename%3DNA12878_S1.vcf.gz.tbi&response-content-type=application%2Foctet-stream&Signature=nS2QjrLPMsvAgJFWn%2BHjmXfTdWw%3D' -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:102.0) Gecko/20100101 Firefox/102.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' -H 'Accept-Language: de,en-US;q=0.7,en;q=0.3' -H 'Accept-Encoding: gzip, deflate, br' -H 'Referer: https://basespace.illumina.com/' -H 'Connection: keep-alive' -H 'Upgrade-Insecure-Requests: 1' -H 'Sec-Fetch-Dest: document' -H 'Sec-Fetch-Mode: navigate' -H 'Sec-Fetch-Site: cross-site' -H 'Sec-Fetch-User: ?1' \
--output NA12878_S1.vcf.gz.tbi
curl 'https://basespace-data-east.s3-external-1.amazonaws.com/47c60b82082847db8cf4a7d900e530b4/NA12878_S1.ROH.bed?AWSAccessKeyId=AKIARPYQJSWQRDJKAWUT&Expires=1665568438&response-content-disposition=attachment%3Bfilename%3DNA12878_S1.ROH.bed&response-content-type=application%2Foctet-stream&Signature=yijtINIcFla%2FLMGZPCKLospTIWI%3D' -H 'User-Agent: Mozilla/5.0 (Macintosh; Intel Mac OS X 10.15; rv:102.0) Gecko/20100101 Firefox/102.0' -H 'Accept: text/html,application/xhtml+xml,application/xml;q=0.9,image/avif,image/webp,*/*;q=0.8' -H 'Accept-Language: de,en-US;q=0.7,en;q=0.3' -H 'Accept-Encoding: gzip, deflate, br' -H 'Referer: https://basespace.illumina.com/' -H 'Connection: keep-alive' -H 'Upgrade-Insecure-Requests: 1' -H 'Sec-Fetch-Dest: document' -H 'Sec-Fetch-Mode: navigate' -H 'Sec-Fetch-Site: cross-site' -H 'Sec-Fetch-User: ?1' \
--output NA12878_S1.ROH.bed

# HG005 - NA24631: (there is no truth set so far)
# ${WGET} ${NCBI}/${ReferenceSamples}/...     .vcf.gz
# ${WGET} ${NCBI}/${ReferenceSamples}/...     .vcf.gz.tbi
# ${WGET} ${NCBI}/${ReferenceSamples}/...     .bed

cd ../..
echo "$(tput setaf 1)$(tput setab 7)------- truth set downloaded (3.4/3.4) --------$(tput sgr 0)" 1>&3
