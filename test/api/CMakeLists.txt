cmake_minimum_required (VERSION 3.8)

add_api_test (junction_detection_test.cpp)
target_use_datasources (junction_detection_test FILES simulated.minimap2.hg19.coordsorted_cutoff.sam)

add_api_test (junction_detection_methods_test.cpp)
target_use_datasources (junction_detection_methods_test FILES simulated.minimap2.hg19.coordsorted_cutoff.sam)
