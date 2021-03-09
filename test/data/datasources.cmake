cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# Get Hash-value via bash: 'shasum -a 256 path/to/file.ending

# copies file to <build>/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
# This file is a small cutoff created by:
# head -20999 simulated.minimap2.hg19.coordsorted.sam | tail -4 > simulated.minimap2.hg19.coordsorted_cutoff.sam
declare_datasource (FILE simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL_HASH SHA256=2423106a8d7d72508ab9009a07f2dc814fa3b01b3983d5975339004fd4b1f9df)

# copies file to <build>/data/detect_breakends_shorted.vcf
declare_datasource (FILE detect_breakends_shorted.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/detect_breakends_shorted.vcf
                    URL_HASH SHA256=c81247543fb3603772ddd77f00e30f74d6f05e9e7f5f008ddcd516e18b040dab)

# copies file to <build>/data/detect_breakends_insertion_file_test.fasta
declare_datasource (FILE detect_breakends_insertion_file_test.fasta
                    URL ${CMAKE_SOURCE_DIR}/test/data/detect_breakends_insertion_file_test.fasta
                    URL_HASH SHA256=8b7ffa31b84a829c19dcd9fb4cabfb49efce91eb2498131199cfaea1daedb0dc)
