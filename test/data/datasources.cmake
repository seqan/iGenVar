cmake_minimum_required (VERSION 3.11)

include (cmake/app_datasources.cmake)

# Get Hash-value via bash: 'shasum -a 256 path/to/file.ending

# copies file to <build>/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
# This file is a small cutoff created by:
# head -20999 simulated.minimap2.hg19.coordsorted.sam | tail -4 > simulated.minimap2.hg19.coordsorted_cutoff.sam
declare_datasource (FILE simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL_HASH SHA256=e59b42c85ed309faf8b3d2f1a8e64a9ccd0a47becd1cb291144efd56be0aa4f9)

# creates the link <build>/data/mini_example_reference.fasta
declare_datasource (FILE mini_example_reference.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/mini_example_reference.fasta
                    URL_HASH SHA256=068f60027a7300e4bb342e34e86ea774c1a0a288d6819d5f75341adcaf6edd9e)

# copies file to <build>/data/paired_end_mini_example.sam
declare_datasource (FILE paired_end_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/paired_end_mini_example.sam
                    URL_HASH SHA256=9dc47068c9a685ae414d9d943b512ab14d1d8de041d805e61c7f8352831a51b3)

# copies file to <build>/data/single_end_mini_example.sam
declare_datasource (FILE single_end_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/single_end_mini_example.sam
                    URL_HASH SHA256=5dbf1d7f41b392bd34ff765b65ac9be73b08246aad774a2399a83523ca45cf41)

# copies file to <build>/data/output_err.txt
declare_datasource (FILE output_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_err.txt
                    URL_HASH SHA256=587d7d853a713cc6cb1af211a95190e6bd50b261160f7222d318c0841e9f16ff)

# copies file to <build>/data/output_res.vcf
declare_datasource (FILE output_res.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_res.vcf
                    URL_HASH SHA256=7dd5c100edaecf8227e92825aeb49649a8b1c3f459d047d67687b36b6f7cdaa9)

# copies file to <build>/data/output_short_and_long_err.txt
declare_datasource (FILE output_short_and_long_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_err.txt
                    URL_HASH SHA256=a4aa69df5a9adfd6573e33686fb2c434cea21775420869c08070e458a5786714)

# copies file to <build>/data/output_short_and_long_res.vcf
declare_datasource (FILE output_short_and_long_res.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_res.vcf
                    URL_HASH SHA256=6ca496cf1c10699cf93c2c32dad2005d3ee5a4c9fd1a391c7eed81abd0b74376)
