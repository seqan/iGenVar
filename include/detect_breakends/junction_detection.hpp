#pragma once

#include <seqan3/io/sequence_file/all.hpp>  // FASTA support
#include <seqan3/io/alignment_file/all.hpp> // SAM/BAM support

#include "detect_breakends/bam_functions.hpp"   // for hasFlag* functions
#include "junction.hpp"                         // for class junction
#include "detect_breakends/aligned_segment.hpp" // for struct aligned_segment

void detect_junctions_in_alignment_file(const std::filesystem::path & alignment_file_path,
                                        const std::filesystem::path & insertion_file_path);
