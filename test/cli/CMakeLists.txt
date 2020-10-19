cmake_minimum_required (VERSION 3.8)

add_cli_test(find_deletions_cli_test.cpp)
add_dependencies (find_deletions_cli_test "find_deletions")
target_use_datasources (find_deletions_cli_test FILES detect_breakends_shorted.vcf)

add_cli_test(detect_breakends_cli_test.cpp)
add_dependencies (detect_breakends_cli_test "detect_breakends")
target_use_datasources (detect_breakends_cli_test FILES simulated.minimap2.hg19.coordsorted_cutoff.sam)
target_use_datasources (detect_breakends_cli_test FILES detect_breakends_insertion_file_test.fasta)

# add_cli_test (iGenVar_options_test.cpp)
# target_use_datasources (iGenVar_options_test FILES in.fastq)
