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
                    URL_HASH SHA256=5e9211c3cb68105bd3359a1f001c900489bd2a5d117fa0d0da82e83766a8325d)

# copies file to <build>/data/single_end_mini_example.sam
declare_datasource (FILE single_end_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/single_end_mini_example.sam
                    URL_HASH SHA256=2c478cade5851a2b8a3aa7db0064be29bfbe137c8844c25848a6b7698ce69f0a)

# copies file to <build>/data/output_err.txt
declare_datasource (FILE output_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_err.txt
                    URL_HASH SHA256=b4c832bbf50cf3b9893191caa4dcc811299ec9de0270f2470bbb22a509826d6a)

# copies file to <build>/data/output_res.txt
declare_datasource (FILE output_res.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_res.txt
                    URL_HASH SHA256=6ce03518859776b951b9246c0bdf62537484cdb120de752468597f21524523de)

# copies file to <build>/data/output_short_and_long_err.txt
declare_datasource (FILE output_short_and_long_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_err.txt
                    URL_HASH SHA256=8dcede8d8b733e1ac9eb2b327c5d14d890d0c6f4e77c95ec8bf9e0d10a593a2e)

# copies file to <build>/data/output_short_and_long_res.txt
declare_datasource (FILE output_short_and_long_res.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_res.txt
                    URL_HASH SHA256=d17b7b163bdbbb2c36fab6ca8a6282645e1ceaf240f66c183c56ccb307d59b90)
