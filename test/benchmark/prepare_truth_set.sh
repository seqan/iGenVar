cd truth_set/

TRUTH_SET=HG002_SVs_Tier1_v0.6

## VCF
# unzip vcf
bgzip -d -c ${TRUTH_SET}.vcf.gz > ${TRUTH_SET}.vcf
# add chr to chromosome names
sed -e 's/##contig=<ID=/##contig=<ID=chr/g' ${TRUTH_SET}.vcf > ${TRUTH_SET}.renamed_contig.vcf
awk '{if($0q!~ /^#/) print "chr"$0; else print $0}' ${TRUTH_SET}.renamed_contig.vcf > ${TRUTH_SET}.renamed_chr.vcf
# zip file again
bgzip -c ${TRUTH_SET}.renamed_chr.vcf > ${TRUTH_SET}.renamed_chr.vcf.gz
# create tbi files
tabix -p vcf ${TRUTH_SET}.vcf.gz
tabix -p vcf ${TRUTH_SET}.renamed_chr.vcf.gz

## BED
# add chr to chromosome names
sed -e 's/^/chr/' ${TRUTH_SET}.bed > ${TRUTH_SET}.renamed_chr.bed

## Create splitted VCF:
mkdir -p INS/
mkdir -p DEL/
# filter for SVTYPES
less ${TRUTH_SET}.vcf | head -102 > INS/${TRUTH_SET}.INS.vcf
less ${TRUTH_SET}.vcf | grep "SVTYPE=INS" >> INS/${TRUTH_SET}.INS.vcf
less ${TRUTH_SET}.vcf | head -102 > DEL/${TRUTH_SET}.DEL.vcf
less ${TRUTH_SET}.vcf | grep "SVTYPE=DEL" >> DEL/${TRUTH_SET}.DEL.vcf
less ${TRUTH_SET}.renamed_chr.vcf | head -102 > INS/${TRUTH_SET}.INS.renamed_chr.vcf
less ${TRUTH_SET}.renamed_chr.vcf | grep "SVTYPE=INS" >> INS/${TRUTH_SET}.INS.renamed_chr.vcf
less ${TRUTH_SET}.renamed_chr.vcf | head -102 > DEL/${TRUTH_SET}.DEL.renamed_chr.vcf
less ${TRUTH_SET}.renamed_chr.vcf | grep "SVTYPE=DEL" >> DEL/${TRUTH_SET}.DEL.renamed_chr.vcf
# zip new files
bgzip -c INS/${TRUTH_SET}.INS.vcf > INS/${TRUTH_SET}.INS.vcf.gz
bgzip -c DEL/${TRUTH_SET}.DEL.vcf > DEL/${TRUTH_SET}.DEL.vcf.gz
bgzip -c INS/${TRUTH_SET}.INS.renamed_chr.vcf > INS/${TRUTH_SET}.INS.renamed_chr.vcf.gz
bgzip -c DEL/${TRUTH_SET}.DEL.renamed_chr.vcf > DEL/${TRUTH_SET}.DEL.renamed_chr.vcf.gz
# index files (create tbi files)
tabix -p vcf INS/${TRUTH_SET}.INS.vcf.gz
tabix -p vcf DEL/${TRUTH_SET}.DEL.vcf.gz
tabix -p vcf INS/${TRUTH_SET}.INS.renamed_chr.vcf.gz
tabix -p vcf DEL/${TRUTH_SET}.DEL.renamed_chr.vcf.gz

echo "$(tput setaf 1)$(tput setab 7)------- truth set files prepared (5.6/5.6) --------$(tput sgr 0)" 1>&3

cd ..
