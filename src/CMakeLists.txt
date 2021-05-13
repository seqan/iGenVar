cmake_minimum_required (VERSION 3.11)

# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib" STATIC modules/clustering/hierarchical_clustering_method.cpp
                                          modules/clustering/simple_clustering_method.cpp
                                          modules/sv_detection_methods/analyze_cigar_method.cpp
                                          modules/sv_detection_methods/analyze_read_pair_method.cpp
                                          modules/sv_detection_methods/analyze_sa_tag_method.cpp
                                          structures/aligned_segment.cpp
                                          structures/breakend.cpp
                                          structures/cluster.cpp
                                          structures/junction.cpp
                                          variant_detection/method_enums.cpp
                                          variant_detection/variant_detection.cpp
                                          variant_detection/variant_output.cpp)

target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC seqan3::seqan3)
target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC fastcluster)
target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)

add_executable ("${PROJECT_NAME}" iGenVar.cpp)
target_link_libraries ("${PROJECT_NAME}" PRIVATE "${PROJECT_NAME}_lib")
