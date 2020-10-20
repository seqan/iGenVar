cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# Get Hash-value via bash: 'shasum -a 256 path/to/file.ending

# copies file to <build>/data/converted_bam_shorted.sam
declare_datasource (FILE converted_bam_shorted.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/converted_bam_shorted.sam
                    URL_HASH SHA256=cde4647e080dda882a0748a0ec34613b8de150a09b90cd831f85034003a9828c)

# copies file to <build>/data/detect_breakends_shorted.vcf
declare_datasource (FILE detect_breakends_shorted.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/detect_breakends_shorted.vcf
                    URL_HASH SHA256=c81247543fb3603772ddd77f00e30f74d6f05e9e7f5f008ddcd516e18b040dab)
