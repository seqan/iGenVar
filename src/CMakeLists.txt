cmake_minimum_required (VERSION 3.8)

# An object library (without main) to be used in multiple targets.
add_library ("${PROJECT_NAME}_lib" STATIC detect_breakends/junction_detection.cpp
                                          find_deletions/deletion_finding_and_printing.cpp
                                          modules/clustering/simple_clustering_method.cpp
                                          modules/sv_detecting_methods/analyze_cigar_method.cpp
                                          modules/sv_detecting_methods/analyze_sa_tag_method.cpp
                                          structures/aligned_segment.cpp
                                          structures/breakend.cpp
                                          structures/cluster.cpp
                                          structures/junction.cpp)

target_link_libraries ("${PROJECT_NAME}_lib" PUBLIC seqan3::seqan3)
target_include_directories ("${PROJECT_NAME}_lib" PUBLIC ../include)

# add_executable ("${PROJECT_NAME}" main.cpp)
# target_link_libraries ("${PROJECT_NAME}" "${PROJECT_NAME}_lib")

add_executable (detect_breakends detect_breakends.cpp)
target_link_libraries (detect_breakends "${PROJECT_NAME}_lib")

add_executable (find_deletions find_deletions.cpp)
target_link_libraries (find_deletions "${PROJECT_NAME}_lib")
