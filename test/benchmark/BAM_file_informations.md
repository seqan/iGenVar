# GRCh37

## Short reads

### HG002.hs37d5.2x250.bam - NIST Illumina 2x250bps Paired-end

#### samtools view -H HG002.hs37d5.2x250.bam

@HD     VN:1.0  SO:unsorted
@PG     ID:novoalign    PN:novoalign    VN:V3.02.07     CL:novoalign -d /cluster/ifs/projects/Genomes/GIAB/refseqs/hs37d5.ndx -f ../../fastq_2x250/D1_S1_L001_R1_001.fastq.gz ../../fastq_2x250/D1_S1_L001_R2_001.fastq.gz -F STDFQ --Q2Off -t 700 -o SAM -c 10
@PG     ID:samtools     PN:samtools     PP:novoalign    VN:1.14 CL:samtools view -H HG002.hs37d5.2x250.bam
@SQ     SN:1    LN:249250621    AS:hs37d5
...
@SQ     SN:22   LN:51304566     AS:hs37d5
@SQ     SN:X    LN:155270560    AS:hs37d5
@SQ     SN:Y    LN:59373566     AS:hs37d5
@SQ     SN:MT   LN:16569        AS:hs37d5
@SQ     SN:GL000207.1   LN:4262 AS:hs37d5
@SQ     SN:GL000226.1   LN:15008        AS:hs37d5
...
@SQ     SN:GL000192.1   LN:547496       AS:hs37d5
@SQ     SN:NC_007605    LN:171823       AS:hs37d5
@SQ     SN:hs37d5       LN:35477943     AS:hs37d5


### HG002.mate_pair.sorted.bam (GRCh37)
from ${NCBI}/AshkenazimTrio/HG002_NA24385_son/NIST_Stanford_Illumina_6kb_matepair/bams/HG002.mate_pair.sorted.bam

#### README
Reads were then mapped to the hg19 reference genome from ucsc
or the GRCh38 reference genome with decoy but no alts using bwa mem (Li 2013) with default settings, and duplicates were marked using samblaster (Faust 2014)
(https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Stanford_Illumina_6kb_matepair/README.NIST_Stanford_Illumina_6kb_matepair)

#### samtools view -H HG002.mate_pair.sorted.bam

@HD     VN:1.3  SO:coordinate
@SQ     SN:chrM LN:16571
@SQ     SN:chr1 LN:249250621
@SQ     SN:chr2 LN:243199373
@SQ     SN:chr3 LN:198022430
@SQ     SN:chr4 LN:191154276
@SQ     SN:chr5 LN:180915260
@SQ     SN:chr6 LN:171115067
@SQ     SN:chr7 LN:159138663
@SQ     SN:chr8 LN:146364022
@SQ     SN:chr9 LN:141213431
@SQ     SN:chr10        LN:135534747
@SQ     SN:chr11        LN:135006516
@SQ     SN:chr12        LN:133851895
@SQ     SN:chr13        LN:115169878
@SQ     SN:chr14        LN:107349540
@SQ     SN:chr15        LN:102531392
@SQ     SN:chr16        LN:90354753
@SQ     SN:chr17        LN:81195210
@SQ     SN:chr18        LN:78077248
@SQ     SN:chr19        LN:59128983
@SQ     SN:chr20        LN:63025520
@SQ     SN:chr21        LN:48129895
@SQ     SN:chr22        LN:51304566
@SQ     SN:chrX LN:155270560
@SQ     SN:chrY LN:59373566
@PG     ID:bwa  PN:bwa  VN:0.7.7-r441   CL:bwa mem -t 2 ref.fa left.fq.gz right.fq.gz
@PG     ID:samtools     PN:samtools     PP:bwa  VN:1.14 CL:samtools view -H HG002.mate_pair.sorted.bam

## Long reads

#### samtools view -H data/long_reads/GRCh37/HG002_PacBio_GRCh37.bam

@HD     VN:1.3  SO:coordinate
@SQ     SN:1    LN:249250621
...
@SQ     SN:22   LN:51304566
@SQ     SN:X    LN:155270560
@SQ     SN:Y    LN:59373566
@SQ     SN:MT   LN:16569
@SQ     SN:GL000207.1   LN:4262
@SQ     SN:GL000226.1   LN:15008
...
@SQ     SN:GL000192.1   LN:547496
@SQ     SN:NC_007605    LN:171823
@SQ     SN:hs37d5       LN:35477943
@RG     ID:HG002-m141008_121115_42156   SM:HG002
...
@RG     ID:HG002-m150124_201043_42177   SM:HG002
@PG     ID:minimap2     PN:minimap2     VN:2.11-r797    CL:minimap2 -t 8 -x map-pb -a --eqx -L -O 5,56 -E 4,1 -B 5 --secondary=no -z 400,50 -r 2k -Y hs37d5.fa -
@PG     ID:samtools     PN:samtools     PP:minimap2     VN:1.14 CL:samtools view -H data/long_reads/GRCh37/HG002_PacBio_GRCh37.bam

#### samtools view -H data/long_reads/GRCh37/HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam

@HD     VN:1.5  SO:coordinate   pb:3.0.1
@SQ     SN:1    LN:249250621    M5:1b22b98cdeb4a9304cb5d48026a85128
@SQ     SN:2    LN:243199373    M5:a0d9851da00400dec1098a9255ac712e
...
@SQ     SN:22   LN:51304566     M5:a718acaa6135fdca8357d5bfe94211dd
@SQ     SN:X    LN:155270560    M5:7e0e2e580297b7764e31dbc80c2540dd
@SQ     SN:Y    LN:59373566     M5:1fa3474750af0948bdf97d5a0ee52e51
@SQ     SN:MT   LN:16569        M5:c68f52674c9fb33aef52dcf399755519
@SQ     SN:GL000207.1   LN:4262 M5:f3814841f1939d3ca19072d9e89f3fd7
@SQ     SN:GL000226.1   LN:15008        M5:1c1b2cd1fccbc0a99b6a447fa24d1504
...
@SQ     SN:GL000192.1   LN:547496       M5:325ba9e808f669dfeee210fdd7b470ac
@SQ     SN:NC_007605    LN:171823       M5:6743bd63b3ff2b5b8985d8933c53290a
@SQ     SN:hs37d5       LN:35477943     M5:5b6a4b3a81a2d3c134b7d14bf6ad39f1
@RG     ID:35690b08     SM:HG002        PU:m54238_180628_014238 DS:READTYPE=CCS;BINDINGKIT=101-365-900;SEQUENCINGKIT=101-309-400;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000 PL:PACBIO       PM:SEQUEL
...
@RG     ID:bbb755dc     SM:HG002        PU:m54316_180719_115850 DS:READTYPE=CCS;BINDINGKIT=101-365-900;SEQUENCINGKIT=101-309-400;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000 PL:PACBIO       PM:SEQUEL
@PG     PN:ccs  ID:ccs-3.0.0    VN:3.0.0        DS:Generate circular consensus sequences (ccs) from subreads.   CL:ccs
...
@PG     PN:ccs  ID:ccs-3.0.0-26EDC745   VN:3.0.0        DS:Generate circular consensus sequences (ccs) from subreads.   CL:ccs  @PG     PN:pbmm2        ID:pbmm2-1C7A8BD3       VN:0.10.0 (commit v0.9.0-25-g3cc3fd9)   CL:pbmm2 align --alignment-threads 16 --preset CCS --sort --sample-name HG002 /pbi/dept/appslab/projects/2018/wr_human_wgs_ccs_kiwi/filteredconsensusreadsets_qv20/m54316_180719_115850.consensusreadset.filtered.xml /pbi/dept/appslab/datasets/wr_references/hs37d5/this.fasta /pbi/dept/appslab/projects/2018/wr_human_wgs_ccs_kiwi/filteredconsensusreadsets_qv20/m54316_180719_115850.Q20.hs37d5.pbmm2.consensusalignmentset.xml
@PG     PN:whatshap     ID:whatshap     VN:0.17 CL:whatshap haplotag --reference /pbi/dept/appslab/datasets/wr_references/hs37d5/this.fasta --output HG002.Sequel.10kb.pbmm2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam ../../RTG.hg19.10x.trio-whatshap.vcf.gz /pbi/dept/appslab/projects/old/2018/wr_hg002_data_sharing/hg19_pbmm2/HG002.Q20.hs37d5.pbmm2.bam  m5:9d47df1972afe81d70dc4e4c59669352
@PG     ID:samtools     PN:samtools     PP:whatshap     VN:1.14 CL:samtools view -H data/long_reads/GRCh37/HG002.Sequel.10kb.pbm
m2.hs37d5.whatshap.haplotag.RTG.10x.trio.bam

#### samtools view -H data/long_reads/GRCh37/NA24385_phased_possorted_bam.bam

@HD     VN:1.3  SO:coordinate
@SQ     SN:chr1 LN:249250621
@SQ     SN:chr2 LN:243199373
...
@SQ     SN:chr22        LN:51304566
@SQ     SN:chrX LN:155270560
@SQ     SN:chrY LN:59373566
@SQ     SN:chrM LN:16571
@SQ     SN:chr1_gl000191_random LN:106433
...
@SQ     SN:chrUn_gl000249       LN:38502
@RG     ID:None SM:10353
@PG     PN:bwa  ID:bwa  VN:0.7.12-r1039 CL:bwa mem -p -t 4 -M -R @RG\tID:None\tSM:10353 /mnt/opt/refdata/fasta/hg19/hg19.fa /mnt/flash2/marsoc/C7A29ANXX.PHASER_SVCALLER_PD.10353.1006.1.0-0/PHASER_SVCALLER_PD/PHASER_SVCALLER/_ALIGNER/TRIM_READS/fork0/chnk0/files/read1.fastq
@PG     PN:10X longranger/attach_bcs    ID:attach_bcs   VN:1006.1.0
@PG     PN:10X longranger/mark_duplicates       ID:mark_duplicates      VN:1006.1.0
@PG     PN:10X longranger/attach_phasing        ID:attach_phasing       VN:1006.1.0
@PG     ID:samtools     PN:samtools     PP:attach_phasing       VN:1.14 CL:samtools view -H data/long_reads/GRCh37/NA24385_phased_possorted_bam.bam

# GRCh38

## Short reads

### HG002.GRCh38.2x250.bam - NIST Illumina 2x250bps Paired-end
#### README
This README describes the library preparation and sequencing performed at NIST for the AJ Trio for the 2x250bp overlapping libraries with nominally 350bp insert size, designed for DISCOVAR assembly.  This folder contains the data for HG002, the AJ son, only.  The folder â€œreadsâ€ contains fastq files from one flow cell, containing ~40-50x coverage.
(https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/README_NIST_Illumina_pairedend_2x250_HG002.txt)

READ DATA:
HG002_NA24385_son/NIST_Illumina_2x250bps/reads
HG003_NA24149_father/NIST_Illumina_2x250bps/reads
HG004_NA24143_mother/NIST_Illumina_2x250bps/reads
REFERENCE FILES:
<!-- "hs37d5" - The reference file available at ftp://ftp.1000genomes.ebi.ac.uk/vol1/ftp/technical/reference/phase2_reference_assembly_sequence/hs37d5.fa.gz -->
"GRCh38" - The reference file available at ftp://ftp.ncbi.nlm.nih.gov/genomes/all/GCA_000001405.15_GRCh38/seqs_for_alignment_pipelines.ucsc_ids/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna.gz
(https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/NIST_Illumina_2x250bps/novoalign_bams/README_update_feb2019)


#### samtools view -H HG002.GRCh38.2x250.bam
@HD     VN:1.0  SO:unsorted
@PG     ID:novoalign    PN:novoalign    VN:V3.02.07     CL:novoalign -d /cluster/ifs/projects/Genomes/GIAB/refseqs/GRCh38_full_plus_hs38d1_analysis_set_minus_alts.ndx -f ../../
fastq_2x250/D1_S1_L001_R1_001.fastq.gz ../../fastq_2x250/D1_S1_L001_R2_001.fastq.gz -F STDFQ --Q2Off -t 700 -o SAM -c 10
@PG     ID:samtools     PN:samtools     PP:novoalign    VN:1.14 CL:samtools view -H HG002.GRCh38.2x250.bam
@SQ     SN:chr1 LN:248956422    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr2 LN:242193529    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr3 LN:198295559    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr4 LN:190214555    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr5 LN:181538259    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr6 LN:170805979    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr7 LN:159345973    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr8 LN:145138636    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr9 LN:138394717    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr10        LN:133797422    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr11        LN:135086622    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr12        LN:133275309    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr13        LN:114364328    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr14        LN:107043718    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr15        LN:101991189    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr16        LN:90338345     AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr17        LN:83257441     AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr18        LN:80373285     AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr19        LN:58617616     AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr20        LN:64444167     AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr21        LN:46709983     AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr22        LN:50818468     AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chrX LN:156040895    AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chrY LN:57227415     AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chrM LN:16569        AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr1_KI270706v1_random       LN:175055       AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chr1_KI270707v1_random       LN:32032        AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
...
@SQ     SN:chr22_KI270739v1_random      LN:73985        AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chrY_KI270740v1_random       LN:37240        AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chrUn_KI270302v1     LN:2274 AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chrUn_KI270304v1     LN:2165 AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
...
@SQ     SN:chrUn_JTFH01001997v1_decoy   LN:2003 AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts
@SQ     SN:chrUn_JTFH01001998v1_decoy   LN:2001 AS:GRCh38_full_plus_hs38d1_analysis_set_minus_alts

### HG002.sorted.bam - NIST Stanford Illumina 6kb matepair (GRCh38)

#### README
Reads were then mapped to the hg19 reference genome from ucsc
or the GRCh38 reference genome with decoy but no alts using bwa mem (Li 2013) with default settings, and duplicates were marked using samblaster (Faust 2014)
(https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002 NA24385_son/NIST_Stanford_Illumina_6kb_matepair/README.NIST_Stanford_Illumina_6kb_matepair)

#### samtools view -H HG002.sorted.bam
@HD     VN:1.3  SO:coordinate
@SQ     SN:chr1 LN:248956422
@SQ     SN:chr2 LN:242193529
@SQ     SN:chr3 LN:198295559
@SQ     SN:chr4 LN:190214555
@SQ     SN:chr5 LN:181538259
@SQ     SN:chr6 LN:170805979
@SQ     SN:chr7 LN:159345973
@SQ     SN:chr8 LN:145138636
@SQ     SN:chr9 LN:138394717
@SQ     SN:chr10        LN:133797422
@SQ     SN:chr11        LN:135086622
@SQ     SN:chr12        LN:133275309
@SQ     SN:chr13        LN:114364328
@SQ     SN:chr14        LN:107043718
@SQ     SN:chr15        LN:101991189
@SQ     SN:chr16        LN:90338345
@SQ     SN:chr17        LN:83257441
@SQ     SN:chr18        LN:80373285
@SQ     SN:chr19        LN:58617616
@SQ     SN:chr20        LN:64444167
@SQ     SN:chr21        LN:46709983
@SQ     SN:chr22        LN:50818468
@SQ     SN:chrX LN:156040895
@SQ     SN:chrY LN:57227415
@SQ     SN:chrM LN:16569
@SQ     SN:chr1_KI270706v1_random       LN:175055
@SQ     SN:chr1_KI270707v1_random       LN:32032
...
@SQ     SN:chrUn_KN707991v1_decoy       LN:2193
@SQ     SN:chrUn_KN707992v1_decoy       LN:1830
@SQ     SN:chrUn_JTFH01000001v1_decoy   LN:25139
@SQ     SN:chrUn_JTFH01000002v1_decoy   LN:18532
...
@SQ     SN:chrUn_JTFH01001997v1_decoy   LN:2003
@SQ     SN:chrUn_JTFH01001998v1_decoy   LN:2001
@PG     ID:bwa  PN:bwa  VN:0.7.13-r1126 CL:bwa mem -t 16 /scratch/users/nspies/data/hg38_giab/GCA_000001405.15_GRCh38_no_alt_plus_hs38d1_analysis_set.fna ../HG002.mate_pair.R1.fastq.gz ../HG002.mate_pair.R2.fastq.gz
@PG     ID:samtools     PN:samtools     PP:bwa  VN:1.14 CL:samtools view -H HG002.sorted.bam

## Long reads

### HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam - PacBio CCS 10kb (minimap2)

#### README
reference: reads were aligned to the GRCh38_no_alt_analysis reference using pbmm2
(https://ftp-trace.ncbi.nlm.nih.gov/ReferenceSamples/giab/data/AshkenazimTrio/HG002_NA24385_son/PacBio_CCS_10kb/GRCh38_no_alt_analysis/README.txt)

#### samtools view -H HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam

@HD     VN:1.5  SO:coordinate   pb:3.0.1
@SQ     SN:chr1 LN:248956422
@SQ     SN:chr2 LN:242193529
@SQ     SN:chr3 LN:198295559
@SQ     SN:chr4 LN:190214555
@SQ     SN:chr5 LN:181538259
@SQ     SN:chr6 LN:170805979
@SQ     SN:chr7 LN:159345973
@SQ     SN:chr8 LN:145138636
@SQ     SN:chr9 LN:138394717
@SQ     SN:chr10        LN:133797422
@SQ     SN:chr11        LN:135086622
@SQ     SN:chr12        LN:133275309
@SQ     SN:chr13        LN:114364328
@SQ     SN:chr14        LN:107043718
@SQ     SN:chr15        LN:101991189
@SQ     SN:chr16        LN:90338345
@SQ     SN:chr17        LN:83257441
@SQ     SN:chr18        LN:80373285
@SQ     SN:chr19        LN:58617616
@SQ     SN:chr20        LN:64444167
@SQ     SN:chr21        LN:46709983
@SQ     SN:chr22        LN:50818468
@SQ     SN:chrX LN:156040895
@SQ     SN:chrY LN:57227415
@SQ     SN:chrM LN:16569
@SQ     SN:chr1_KI270706v1_random       LN:175055
@SQ     SN:chr1_KI270707v1_random       LN:32032
...
@SQ     SN:chr22_KI270739v1_random      LN:73985
@SQ     SN:chrY_KI270740v1_random       LN:37240
@SQ     SN:chrUn_KI270302v1     LN:2274
@SQ     SN:chrUn_KI270304v1     LN:2165
...
@SQ     SN:chrUn_GL000216v2     LN:176608
@SQ     SN:chrUn_GL000218v1     LN:161147@SQ     SN:chrEBV       LN:171823
@RG     ID:379a27f3     SM:HG002        PU:m54315_180703_170746 DS:READTYPE=CCS;BINDINGKIT=101-365-900;SEQUENCINGKIT=101-309-400;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000P
L:PACBIO        PM:SEQUEL
@RG     ID:697c9e3c     SM:HG002        PU:m54315_180705_094736 DS:READTYPE=CCS;BINDINGKIT=101-365-900;SEQUENCINGKIT=101-309-400;BASECALLERVERSION=5.0.0;FRAMERATEHZ=100.000000PL:PACBIO        PM:SEQUEL
...
@PG     PN:pbmm2        ID:pbmm2-1A9717AC       VN:1.1.0 (commit v1.0.0-38-g46b3940)    CL:pbmm2 align --alignment-threads 14 --sort-threads 2 --sort-memory 16G --preset CCS -L 0.1 -c 0 --sort --sample HG002 /pbi/dept/appslab/projects/2018/wr_human_wgs_ccs_kiwi/filteredconsensusreadsets_qv20/m54315_180704_132638.consensusreadset.filtered.xml /pbi/dept/secondary/siv/references/human_GRCh38_no_alt_analysis_set/sequence/human_GRCh38_no_alt_analysis_set.fasta 9.consensusalignmentset.xml
@PG     PN:whatshap     ID:whatshap-30CBFBAF    VN:0.17 CL:whatshap haplotag --output 9.haplotagged.bam --reference /pbi/dept/secondary/siv/references/human_GRCh38_no_alt_analysis_set/sequence/human_GRCh38_no_alt_analysis_set.fasta ../NA24385.GRCh38.phased_variants.reheadered.vcf.gz 9.bam       m5:277415e7d359807063a0b6dc7dfd6bed
@PG     ID:samtools     PN:samtools     PP:whatshap-30CBFBAF    VN:1.14 CL:samtools view -H HG002.10kb.Sequel.pbmm2.GRCh38.whatshap.haplotag.RTG.10x.trio.bam
