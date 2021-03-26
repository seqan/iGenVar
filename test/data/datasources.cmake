cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# Get Hash-value via bash: 'shasum -a 256 path/to/file.ending

# copies file to <build>/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
# This file is a small cutoff created by:
# head -20999 simulated.minimap2.hg19.coordsorted.sam | tail -4 > simulated.minimap2.hg19.coordsorted_cutoff.sam
declare_datasource (FILE simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL_HASH SHA256=e59b42c85ed309faf8b3d2f1a8e64a9ccd0a47becd1cb291144efd56be0aa4f9)

# copies file to <build>/data/detect_breakends_shorted.vcf
declare_datasource (FILE detect_breakends_shorted.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/detect_breakends_shorted.vcf
                    URL_HASH SHA256=c81247543fb3603772ddd77f00e30f74d6f05e9e7f5f008ddcd516e18b040dab)

# copies file to <build>/data/detect_breakends_insertion_file_test.fasta
declare_datasource (FILE detect_breakends_insertion_file_test.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/detect_breakends_insertion_file_test.fasta
                    URL_HASH SHA256=8b7ffa31b84a829c19dcd9fb4cabfb49efce91eb2498131199cfaea1daedb0dc)

# copies file to <build>/data/mini_example.sam
declare_datasource (FILE mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/mini_example.sam
                    URL_HASH SHA256=8a8878c7f067770f8d64242c9c66f1f4206801899d843205fa4d5d4da455d0e7)

# copies file to <build>/data/output_err.txt
declare_datasource (FILE output_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_err.txt
                    URL_HASH SHA256=376894af7c15e72a1636a08e593f9c0db87dab7fa142e9c5bfb9be41574e9033)

# copies file to <build>/data/insertions_output.fasta
declare_datasource (FILE insertions_output.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/insertions_output.fasta
                    URL_HASH SHA256=ae47d421f6f7ee33ca9e1d1cb072182b8f39319fbab9802866a5df8e6b281d27)
