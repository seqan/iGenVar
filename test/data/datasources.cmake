cmake_minimum_required (VERSION 3.11)

include (cmake/app_datasources.cmake)

# Get Hash-value via bash: 'shasum -a 256 path/to/file.ending

# copies file to <build>/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
# This file is a small cutoff created by:
# head -20999 simulated.minimap2.hg19.coordsorted.sam | tail -4 > simulated.minimap2.hg19.coordsorted_cutoff.sam
declare_datasource (FILE simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL_HASH SHA256=e59b42c85ed309faf8b3d2f1a8e64a9ccd0a47becd1cb291144efd56be0aa4f9)

# copies file to <build>/data/paired_end_mini_example.sam
declare_datasource (FILE paired_end_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/paired_end_mini_example.sam
                    URL_HASH SHA256=a4134ef84c4a2d6aa57f00966a69ca3c54036dc3d17f418427ce4fc6bb7e87ff)

# copies file to <build>/data/single_end_mini_example.sam
declare_datasource (FILE single_end_mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/single_end_mini_example.sam
                    URL_HASH SHA256=2c478cade5851a2b8a3aa7db0064be29bfbe137c8844c25848a6b7698ce69f0a)

# copies file to <build>/data/output_err.txt
declare_datasource (FILE output_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_err.txt
                    URL_HASH SHA256=386fc7508d96e8bb55333b5d5e36e749e072bad5dc94b3bb1c90b241910d1b60)

# copies file to <build>/data/output_res.txt
declare_datasource (FILE output_res.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_res.txt
                    URL_HASH SHA256=8508f442391cc7118822641a0dcb3cbed3da192482f5cbd1457bb28e03097ddd)
