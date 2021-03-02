#pragma once

#include <seqan3/std/filesystem>    // for filesystem
#include <vector>

#include "method_enums.hpp"         // for enum clustering_methods and refinement_methods

/*! \brief Detects junctions between distant genomic positions by analyzing an alignment file (sam/bam). The detected
 *         junctions are printed on stdout and insertion alleles are stored in a fasta file.
 *
 * \param alignment_file_path input file - path to the sam/bam file
 * \param insertion_file_path output file - path for the fasta file
 * \param methods - list of methods for detecting junctions (0: cigar_string,
 *                                                           1: split_read,
 *                                                           2: read_pairs,
 *                                                           3: read_depth)
 * \param clustering_method method for clustering junctions (0: simple_clustering
 *                                                           1: hierarchical_clustering,
 *                                                           2: self-balancing_binary_tree,
 *                                                           3: candidate_selection_based_on_voting)
 * \param refinement_method method for refining breakends (0: no_refinement,
 *                                                         1: sViper_refinement_method,
 *                                                         2: sVirl_refinement_method)
 * \param min_var_length - minimum length of variants to detect (default 30 bp)
 *
 * \details Detects junctions from the CIGAR strings and supplementary alignment tags of read alignment records.
 *          We sort out unmapped alignments, secondary alignments, duplicates and alignments with low mapping quality.
 *          Then, the CIGAR string of all remaining alignments is analyzed.
 *          For primary alignments, also the split read information is analyzed.
 */
void detect_junctions_in_alignment_file(const std::filesystem::path & alignment_file_path,
                                        const std::filesystem::path & insertion_file_path,
                                        const std::vector<detecting_methods> methods,
                                        const clustering_methods clustering_method,
                                        const refinement_methods refinement_method,
                                        const uint64_t min_var_length);
