cmake_minimum_required (VERSION 3.16)

include (cmake/app_datasources.cmake)

# Get Hash-value via bash: 'shasum -a 256 path/to/file.ending

# copies file to <build>/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
# This file is a small cutoff created by:
# head -20999 simulated.minimap2.hg19.coordsorted.sam | tail -4 > simulated.minimap2.hg19.coordsorted_cutoff.sam
declare_datasource (FILE simulated.minimap2.hg19.coordsorted_cutoff.sam URL
                    ${CMAKE_SOURCE_DIR}/test/data/simulated.minimap2.hg19.coordsorted_cutoff.sam URL_HASH
                    SHA256=e59b42c85ed309faf8b3d2f1a8e64a9ccd0a47becd1cb291144efd56be0aa4f9)

# creates the link <build>/data/mini_example_reference.fasta
declare_datasource (FILE mini_example_reference.fasta URL
                    ${CMAKE_SOURCE_DIR}/test/data/mini_example/mini_example_reference.fasta URL_HASH
                    SHA256=068f60027a7300e4bb342e34e86ea774c1a0a288d6819d5f75341adcaf6edd9e)

# copies file to <build>/data/paired_end_mini_example.sam
declare_datasource (FILE paired_end_mini_example.sam URL
                    ${CMAKE_SOURCE_DIR}/test/data/mini_example/paired_end_mini_example.sam URL_HASH
                    SHA256=58961c8f016dbcaa6d38806f9cd3e90dbacda4f40af4693db9540c5544b9367a)

# copies file to <build>/data/single_end_mini_example.sam
declare_datasource (FILE single_end_mini_example.sam URL
                    ${CMAKE_SOURCE_DIR}/test/data/mini_example/single_end_mini_example.sam URL_HASH
                    SHA256=b7ad2c43c444e4897e883c01f8efe3bbdd2b185d1c30efaa16be28ca2b04049d)

# copies file to <build>/data/output_err.txt
declare_datasource (FILE output_err.txt URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_err.txt URL_HASH
                    SHA256=d30107b5222a0768db9cbe6f4f7c4b11b7c530057f99c969d891a56f596f6320)

# copies file to <build>/data/output_res.vcf
declare_datasource (FILE output_res.vcf URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_res.vcf URL_HASH
                    SHA256=938d73ef495f4cddccc10eaefbda7c3ec934926798edf5a8db4d1b00ae10e646)

# copies file to <build>/data/output_short_and_long_err.txt
declare_datasource (FILE output_short_and_long_err.txt URL
                    ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_err.txt URL_HASH
                    SHA256=86dfe165ea871101b52c37dbcb6db8b4ccf86a2484d12277ecdd4b5e75be5cc8)

# copies file to <build>/data/output_short_and_long_res.vcf
declare_datasource (FILE output_short_and_long_res.vcf URL
                    ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_res.vcf URL_HASH
                    SHA256=a95fd7e6f882ea977feaef541ec8fc92906b5fda329088433670fcaa05d36a50)
