cmake_minimum_required (VERSION 3.11)

add_cli_test(iGenVar_cli_test.cpp)
add_dependencies (iGenVar_cli_test "iGenVar")
target_use_datasources (iGenVar_cli_test FILES simulated.minimap2.hg19.coordsorted_cutoff.sam)
target_use_datasources (iGenVar_cli_test FILES mini_example.sam)
target_use_datasources (iGenVar_cli_test FILES output_err.txt)

# add_cli_test (iGenVar_options_test.cpp)
# target_use_datasources (iGenVar_options_test FILES in.fastq)
