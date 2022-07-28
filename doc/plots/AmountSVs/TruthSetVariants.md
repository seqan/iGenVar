`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=" | wc -l`
74012
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=" | wc -l`
74012
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=" | grep "REPTYPE=" | wc -l`
74012

# INS
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INS" | wc -l`
36600
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INS" | grep "REPTYPE=SIMPLEINS" | wc -l`
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=SIMPLEINS" | wc -l`
7639
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=DUP" | wc -l`
0
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=DUP" | wc -l`
27309
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INS" | grep "REPTYPE=DUP" | wc -l`
27297
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INS" | grep "REPTYPE=SIMPLEDEL" | wc -l`
234
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=SUBSINS" | wc -l`
1263
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INS" | grep "REPTYPE=SUBSINS" | wc -l`
1261
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INS" | grep "REPTYPE=SUBSDEL" | wc -l`
78
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INS" | grep "REPTYPE=CONTRAC" | wc -l`
90

# DEL
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=DEL" | wc -l`
37412
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=SIMPLEDEL" | wc -l`
24765
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=DEL" | grep "REPTYPE=SIMPLEDEL" | wc -l`
24531
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=SUBSDEL" | wc -l`
1077
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=DEL" | grep "REPTYPE=SUBSINS" | wc -l`
2
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=DEL" | grep "REPTYPE=SUBSDEL" | wc -l`
999
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=DEL" | grep "REPTYPE=DUP" | wc -l`
12
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=CONTRAC" | wc -l`
11957
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=DEL" | grep "REPTYPE=CONTRAC" | wc -l`
11867

`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=INV" | wc -l`
0
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=CON" | wc -l`
0
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INV" | wc -l`
0

# SVLEN
## SV < 20bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1<20{print $1}' | wc -l`
29
## 20bp < SV < 49bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>20&&$1<49{print $1}' | wc -l`
15833
## 50bp < SV < 249bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>50&&$1<249{print $1}' | wc -l`
6829
## 250bp < SV < 499bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>250&&$1<499{print $1}' | wc -l`
2664
## 500bp < SV < 999bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>500&&$1<999{print $1}' | wc -l`
1270
## SV > 1000bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>1000{print $1}' | wc -l`
1321




## SV < 10bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1<10{print $1}' | wc -l`
15
## SV < 30bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1<30{print $1}' | wc -l`
10093

## SV < 50bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1<50{print $1}' | wc -l`
16026

## SV < 500bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1<500{print $1}' | wc -l`
25530

## SV > 1000bp
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>1000{print $1}' | wc -l`
1321
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=INS" | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>1000{print $1}' | wc -l`
1321
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "SVTYPE=DEL" | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>1000{print $1}' | wc -l`
0
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=SIMPLEINS" | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>1000{print $1}' | wc -l`
899
`$ less data/truth_set/HG002_SVs_Tier1_v0.6.vcf | grep "REPTYPE=DUP" | grep 'SVLEN=[1-9]*\;' | awk -F 'SVLEN=' '{print $2}' | awk -F ';' '$1>1000{print $1}' | wc -l`
50
