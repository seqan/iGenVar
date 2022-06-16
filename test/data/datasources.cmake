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
                    URL_HASH SHA256=e19a354491cebf80abafb4fbf6886419dab9829bcd55f27adc0ca67ddfe51a1a)

# copies file to <build>/data/output_res.vcf
declare_datasource (FILE output_res.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_res.vcf
                    URL_HASH SHA256=245974501fe0b8cfc63873008780c46eec115ecd63f90809b6fa98f2eb6852ca)

# copies file to <build>/data/output_short_and_long_err.txt
declare_datasource (FILE output_short_and_long_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_err.txt
                    URL_HASH SHA256=1f729eeda96a0cba213a4bec8de9ef6ea1034c1079f645782400bee4e12aced2)

# copies file to <build>/data/output_short_and_long_res.vcf
declare_datasource (FILE output_short_and_long_res.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_res.vcf
                    URL_HASH SHA256=b32b53637e93c769f21fd92345cc5194788f44abd8f021778c2624e917e87b3a)
