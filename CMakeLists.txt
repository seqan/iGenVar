cmake_minimum_required (VERSION 3.22)

## CUSTOMISE

# Define the application name and version.
project (iGenVar VERSION 1.0.0)

## BUILD

# Make Release default build type
if (NOT CMAKE_BUILD_TYPE)
    set (CMAKE_BUILD_TYPE Release CACHE STRING
         "Choose the type of build, options are: Debug Release Coverage RelWithDebInfo MinSizeRel."
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

list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/lib/seqan3/build_system")
list (APPEND CMAKE_PREFIX_PATH "${CMAKE_CURRENT_LIST_DIR}/lib/b.i.o./build_system")

# Dependency: SeqAn3 is found by b.i.o. as well.
find_package (SeqAn3 REQUIRED)

# Dependency: b.i.o.
set (BIO_DEPENDENCY_INCLUDE_DIRS "${CMAKE_CURRENT_LIST_DIR}/lib/seqan3/include")
set (BIO_SUBMODULES_DIR "/dev/null") # Prevent b.i.o. from including its seqan3 submodule
find_package (bio REQUIRED)


# Use ccache.
include ("${SEQAN3_CLONE_DIR}/test/cmake/seqan3_require_ccache.cmake")
seqan3_require_ccache ()

# Dependency: hclust.
include (FetchContent)
FetchContent_Declare(
        hclust
        GIT_REPOSITORY https://github.com/cdalitz/hclust-cpp.git
        GIT_TAG 8d2fe71fdc2dc5e95e76485f54b9f075a63a54e6        # commit from 17 Aug 2021
)
FetchContent_Populate(hclust)
add_library ("fastcluster" STATIC ${hclust_SOURCE_DIR}/fastcluster.cpp)
target_include_directories ("fastcluster" PUBLIC ${hclust_SOURCE_DIR})

# Dependency: Lemon Graph Library
FetchContent_Declare(
        lemon
        GIT_REPOSITORY https://github.com/seqan/lemon.git
        GIT_TAG 400d51cff5b4c4b59dd11a30965b9837a3df7b53
)
FetchContent_GetProperties(lemon)
if (NOT lemon_POPULATED)
    FetchContent_Populate(lemon)
    include_directories(SYSTEM ${lemon_SOURCE_DIR}/include)
endif ()

# Dependency: BAMIntervalTree.
add_library ("bamit" INTERFACE)
target_include_directories ("bamit" INTERFACE lib/BAMIntervalTree/include/)

# Add the application.
add_subdirectory (src)
message (STATUS "${FontBold}You can run `make` to build the application.${FontReset}")

## DOCUMENTATION

add_subdirectory (doc EXCLUDE_FROM_ALL)

## TEST

enable_testing ()
add_subdirectory (test EXCLUDE_FROM_ALL)
