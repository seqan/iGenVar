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
                    URL_HASH SHA256=a45fed19a4ef58bcd1ac2d48e4a30d82dd74d7cadc64e468da2487946ce41afc)

# copies file to <build>/data/paired_end_mini_example.sam
declare_datasource (FILE paired_end_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/paired_end_mini_example.sam
                    URL_HASH SHA256=9dc47068c9a685ae414d9d943b512ab14d1d8de041d805e61c7f8352831a51b3)

# copies file to <build>/data/single_end_mini_example.sam
declare_datasource (FILE single_end_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/single_end_mini_example.sam
                    URL_HASH SHA256=2c478cade5851a2b8a3aa7db0064be29bfbe137c8844c25848a6b7698ce69f0a)

# copies file to <build>/data/output_err.txt
declare_datasource (FILE output_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_err.txt
                    URL_HASH SHA256=e70dc7eeb4008d0a36e0515db4bdc83f3e49cb53295d1127bbd0c463eec77696)

# copies file to <build>/data/output_res.vcf
declare_datasource (FILE output_res.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_res.vcf
                    URL_HASH SHA256=98238d077a44d1dfbe8ee342cdeb3e793ef0f026e54fc3f8bf4fe3d9c29a637d)

# copies file to <build>/data/output_short_and_long_err.txt
declare_datasource (FILE output_short_and_long_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_err.txt
                    URL_HASH SHA256=819b50d57e4c0238f3216e6b982072749b568dcc300e8761c9cd4251b6b9daef)

# copies file to <build>/data/output_short_and_long_res.vcf
declare_datasource (FILE output_short_and_long_res.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_res.vcf
                    URL_HASH SHA256=ac250f321cb893d3d888ad34c94cd487913bf509b5808719bcbd4b88ddb1641d)
