#pragma once

#include <seqan3/std/filesystem>    // for filesystem
#include <vector>

#include "method_enums.hpp"         // for enum clustering_methods and refinement_methods

/*! \brief Detects genomic variants by analyzing an alignment file (sam/bam). The detected
 *         variants are printed to a given file or stdout and insertion alleles are stored in a fasta file.
 *
 * \param alignment_short_reads_file_path - short reads input file, path to the sam/bam file
 * \param alignment_long_reads_file_path - long reads input file, path to the sam/bam file
 * \param insertion_file_path output file - path for the fasta file
 * \param methods - list of methods for detecting junctions (1: cigar_string,
 *                                                           2: split_read,
 *                                                           3: read_pairs,
 *                                                           4: read_depth)
 * \param clustering_method method for clustering junctions (0: simple_clustering
 *                                                           1: hierarchical_clustering,
 *                                                           2: self-balancing_binary_tree,
 *                                                           3: candidate_selection_based_on_voting)
 * \param refinement_method method for refining junctions (0: no_refinement,
 *                                                         1: sViper_refinement_method,
 *                                                         2: sVirl_refinement_method)
 * \param min_var_length - minimum length of variants to detect (default 30 bp)
 * \param output_file_path output file - path for the VCF file
 *
 * \details Detects novel junctions from read alignment records using different detection methods.
 *          The junctions are clustered using one of several clustering methods.
 *          Then, the junction clusters are refined using one of several refinement methods.
 *          Finally, the refined junction clusters are categorized into different variant classes
 *          and output in VCF format.
 */
void detect_variants_in_alignment_file(const std::filesystem::path & alignment_short_reads_file_path,
                                       const std::filesystem::path & alignment_long_reads_file_path,
                                       const std::filesystem::path & insertion_file_path,
                                       const std::vector<detection_methods> & methods,
                                       const clustering_methods & clustering_method,
                                       const refinement_methods & refinement_method,
                                       const uint64_t & min_var_length,
                                       const std::filesystem::path & output_file_path);
