cmake_minimum_required (VERSION 3.8)

add_library ("fastcluster" STATIC hclust/fastcluster.cpp)
target_include_directories ("fastcluster" PUBLIC hclust)