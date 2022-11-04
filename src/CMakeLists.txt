cmake_minimum_required (VERSION 3.16)

# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib" STATIC
             modules/clustering/hierarchical_clustering_method.cpp
             modules/clustering/simple_clustering_method.cpp
             modules/sv_detection_methods/analyze_cigar_method.cpp
             modules/sv_detection_methods/analyze_read_pair_method.cpp
             modules/sv_detection_methods/analyze_split_read_method.cpp
             structures/aligned_segment.cpp
             structures/breakend.cpp
             structures/cluster.cpp
             structures/debruijn_graph.cpp
             structures/junction.cpp
             variant_detection/method_enums.cpp
             variant_detection/snp_indel_detection.cpp
             variant_detection/variant_detection.cpp
             variant_detection/variant_output.cpp)

target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC seqan3::seqan3)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC fastcluster)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC bamit)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC biocpp::core)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC biocpp::io)
target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)
target_compile_options ("${PROJECT_NAME}_lib" PUBLIC "-pedantic" "-Wall" "-Wextra")

add_executable ("${PROJECT_NAME}" iGenVar.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")

# MasonVcfConverter - for Simulation Benchmarks
add_executable (MasonVcfConverter ../test/benchmark/simulation/mason_vcf_converter.cpp)
target_link_libraries (MasonVcfConverter PRIVATE biocpp::core biocpp::io seqan3::seqan3)
