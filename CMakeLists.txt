cmake_minimum_required (VERSION 3.11)

## CUSTOMISE

# Define the application name and version.
project (iGenVar VERSION 1.0.0)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release RelWithDebInfo"
         FORCE)
endif ()

# Specify the directories where to store the built archives, libraries and executables
set (CMAKE_ARCHIVE_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_LIBRARY_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/lib")
set (CMAKE_RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/bin")

# Messages
string (ASCII 27 Esc)
set (FontBold "${Esc}[1m")
set (FontReset "${Esc}[m")

# Dependency: SeqAn3.
find_package (SeqAn3 QUIET REQUIRED HINTS lib/seqan3/build_system)

include ("${SEQAN3_CLONE_DIR}/test/cmake/seqan3_require_ccache.cmake")
seqan3_require_ccache ()

# Dependency: hclust.
include (FetchContent)
FetchContent_Declare(
        hclust
        GIT_REPOSITORY https://github.com/cdalitz/hclust-cpp.git
        GIT_TAG f902638a3ba1be40ffb6ff6a8ddafba91a150c76
)
FetchContent_Populate(hclust)
add_library ("fastcluster" STATIC ${hclust_SOURCE_DIR}/fastcluster.cpp)
target_include_directories ("fastcluster" PUBLIC ${hclust_SOURCE_DIR})

# Code Coverage Configuration
add_library (igenvar_coverage_config INTERFACE)

option (IGENVAR_CODE_COVERAGE "Enable coverage reporting" OFF)
if (IGENVAR_CODE_COVERAGE AND CMAKE_CXX_COMPILER_ID MATCHES "GNU")
  # Add required flags (GCC & LLVM/Clang)
  target_compile_options (igenvar_coverage_config INTERFACE
    --coverage # sets all required flags
  )
  if (CMAKE_BUILD_TYPE AND NOT CMAKE_BUILD_TYPE MATCHES "Debug")
    message (WARNING "Please build code coverage in Debug, i.e., cmake -DCMAKE_BUILD_TYPE=Debug")
  endif ()
  if (CMAKE_VERSION VERSION_GREATER_EQUAL 3.13)
    target_link_options (igenvar_coverage_config INTERFACE --coverage)
  else ()
    target_link_libraries (igenvar_coverage_config INTERFACE --coverage)
  endif ()
endif ()

add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
