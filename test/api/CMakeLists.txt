cmake_minimum_required (VERSION 3.8)

add_api_test (junction_detection_test.cpp)
target_use_datasources (junction_detection_test FILES converted_bam_shorted.sam)

add_api_test (deletion_finding_and_printing_test.cpp)
target_use_datasources (deletion_finding_and_printing_test FILES detect_breakends_shorted.vcf)
