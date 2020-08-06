#pragma once

#include <seqan3/io/sequence_file/all.hpp>  // FASTA support
#include <seqan3/io/alignment_file/all.hpp> // SAM/BAM support

#include "detect_breakends/bam_functions.hpp"
#include "junction.hpp"
#include "detect_breakends/aligned_segment.hpp"

/*! \brief Function, detecting junctions in alignment files (sam/bam).
 *  \param alignment_file_path input file - path to the sam/bam file
 *  \param insertion_file_path output file - path for the fasta file
 */
void detect_junctions_in_alignment_file(const std::filesystem::path & alignment_file_path,
                                        const std::filesystem::path & insertion_file_path);
