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
                    URL_HASH SHA256=f5e4707a89e6a089e87b401aa79547a97f5ecdd1599321409f379d1760e8cf6c)

# copies file to <build>/data/output_res.vcf
declare_datasource (FILE output_res.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_res.vcf
                    URL_HASH SHA256=71238d79330c9efcf61467c20cd998517bf755d64742f60cbbfe0a5399a6c063)

# copies file to <build>/data/output_short_and_long_err.txt
declare_datasource (FILE output_short_and_long_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_err.txt
                    URL_HASH SHA256=cf41b4fac14a870b0be40bbee82e2aa4587148498261821d6eac35e7d532a487)

# copies file to <build>/data/output_short_and_long_res.vcf
declare_datasource (FILE output_short_and_long_res.vcf
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_short_and_long_res.vcf
                    URL_HASH SHA256=55262c77c479cd4d767a43a03d4afd1528d722097de8a81d9e52805b37671e01)
