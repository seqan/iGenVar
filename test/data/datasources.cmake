cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# Get Hash-value via bash: 'shasum -a 256 path/to/file.ending

# copies file to <build>/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
# This file is a small cutoff created by:
# head -20999 simulated.minimap2.hg19.coordsorted.sam | tail -4 > simulated.minimap2.hg19.coordsorted_cutoff.sam
declare_datasource (FILE simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL_HASH SHA256=e59b42c85ed309faf8b3d2f1a8e64a9ccd0a47becd1cb291144efd56be0aa4f9)

# copies file to <build>/data/paired_end_short_read_mini_example.sam
declare_datasource (FILE paired_end_short_read_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/paired_end_short_read_mini_example.sam
                    URL_HASH SHA256=11403c69b9967b345dae78c52c25259f529bae20593db487bffffb54c1efdaff)

# copies file to <build>/data/single_end_mini_example.sam
declare_datasource (FILE single_end_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/single_end_mini_example.sam
                    URL_HASH SHA256=ea89fb00c802e5136ca8dd32362648c8192652a34078bfa25b8ae78695a27cd6)

# copies file to <build>/data/output_err.txt
declare_datasource (FILE output_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_err.txt
                    URL_HASH SHA256=b318fde181d620e8c3ab5e8c9b51b55423eeb578dcd192790654a92109d567cb)
