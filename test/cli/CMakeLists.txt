cmake_minimum_required (VERSION 3.8)

add_cli_test (fastq_to_fasta_options_test.cpp)
target_use_datasources (fastq_to_fasta_options_test FILES in.fastq)

add_cli_test(find_deletions_cli_test.cpp)
add_dependencies (find_deletions_cli_test "find_deletions")
target_use_datasources (find_deletions_cli_test FILES detect_breakends_shorted.vcf)

# add_cli_test (iGenVar_options_test.cpp)
# target_use_datasources (iGenVar_options_test FILES in.fastq)
