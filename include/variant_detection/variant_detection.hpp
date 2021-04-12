#pragma once

#include <seqan3/std/filesystem>    // for filesystem
#include <vector>

#include "method_enums.hpp"         // for enum detection_methods, clustering_methods and refinement_methods
#include "structures/junction.hpp"  // for class Junction

//TODO (irallia): move this back to the iGenVar.cpp with the PR #104
struct cmd_arguments
{
    std::filesystem::path alignment_short_reads_file_path{""};
    std::filesystem::path alignment_long_reads_file_path{""};
    std::filesystem::path output_file_path{};
    std::vector<detection_methods> methods{cigar_string, split_read, read_pairs, read_depth};   // default: all methods
    clustering_methods clustering_method{simple_clustering};                                    // default: simple clustering method
    refinement_methods refinement_method{no_refinement};                                        // default: no refinement
    uint64_t min_var_length = 30;
};

/*! \brief Detects junctions between distant genomic positions by analyzing a long read alignment file (sam/bam). The
 *         detected junctions are stored in a vector.
 *
 * \param[in, out]  junctions - a vector of junctions
 * \param[in]       alignment_long_reads_file_path - long reads input file, path to the sam/bam file
 * \param[in]       methods - list of methods for detecting junctions (0: cigar_string,
 *                                                                     1: split_read,
 *                                                                     2: read_pairs,
 *                                                                     3: read_depth)
 * \param[in]       clustering_method - method for clustering junctions (0: simple_clustering
 *                                                                       1: hierarchical_clustering,
 *                                                                       2: self-balancing_binary_tree,
 *                                                                       3: candidate_selection_based_on_voting)
 * \param[in]       refinement_method - method for refining breakends (0: no_refinement,
 *                                                                     1: sViper_refinement_method,
 *                                                                     2: sVirl_refinement_method)
 * \param[in]       min_var_length - minimum length of variants to detect (default 30 bp)
 *
 * \details Detects junctions from the CIGAR strings and supplementary alignment tags of read alignment records.
 *          We filter unmapped alignments, secondary alignments, duplicates and alignments with low mapping quality.
 *          Then, the CIGAR string of all remaining alignments is analyzed.
 *          For primary alignments, also the split read information is analyzed.
 */
void detect_junctions_in_long_reads_sam_file(std::vector<Junction> & junctions,
                                             std::filesystem::path const & alignment_long_reads_file_path,
                                             std::vector<detection_methods> const & methods,
                                             clustering_methods const & clustering_method,
                                             refinement_methods const & refinement_method,
                                             uint64_t const min_var_length);
