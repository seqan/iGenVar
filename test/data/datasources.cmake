cmake_minimum_required (VERSION 3.8)

include (cmake/app_datasources.cmake)

# Get Hash-value via bash: 'shasum -a 256 path/to/file.ending

# copies file to <build>/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
# This file is a small cutoff created by:
# head -20999 simulated.minimap2.hg19.coordsorted.sam | tail -4 > simulated.minimap2.hg19.coordsorted_cutoff.sam
declare_datasource (FILE simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/simulated.minimap2.hg19.coordsorted_cutoff.sam
                    URL_HASH SHA256=e59b42c85ed309faf8b3d2f1a8e64a9ccd0a47becd1cb291144efd56be0aa4f9)

# copies file to <build>/data/mini_example.sam
declare_datasource (FILE mini_example.sam
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/mini_example.sam
                    URL_HASH SHA256=8a8878c7f067770f8d64242c9c66f1f4206801899d843205fa4d5d4da455d0e7)

# copies file to <build>/data/output_err.txt
declare_datasource (FILE output_err.txt
                    URL ${CMAKE_SOURCE_DIR}/test/data/mini_example/output_err.txt
                    URL_HASH SHA256=10ff325f4997b8dee4bb845df92a5949b770f1574ad8022ef7761d1401f45ed8)
