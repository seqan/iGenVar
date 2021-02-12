#pragma once

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/all.hpp>  // FASTA support
#include <seqan3/io/alignment_file/all.hpp> // SAM/BAM support

#include "detect_breakends/bam_functions.hpp"   // for hasFlag* functions
#include "junction.hpp"                         // for class junction
#include "cluster.hpp"                         // for class cluster
#include "detect_breakends/aligned_segment.hpp" // for struct aligned_segment
#include "method_enums.hpp"                     // for enum clustering_methods

/*! \brief Detects junctions between distant genomic positions by analyzing an alignment file (sam/bam). The detected
 *         junctions are printed on stdout and insertion alleles are stored in a fasta file.
 *
 * \param alignment_file_path input file - path to the sam/bam file
 * \param insertion_file_path output file - path for the fasta file
 * \param methods - list of methods for detecting junctions (1: cigar_string,
 *                                                           2: split_read,
 *                                                           3: read_pairs,
 *                                                           4: read_depth)
 * \param clustering_method method for clustering junctions (0: simple_clustering
 *                                                           1: hierarchical_clustering,
 *                                                           2: self-balancing_binary_tree,
 *                                                           3: candidate_selection_based_on_voting)
 *
 * \details Detects junctions from the CIGAR strings and supplementary alignment tags of read alignment records.
 *          We sort out unmapped alignments, secondary alignments, duplicates and alignments with low mapping quality.
 *          Then, the CIGAR string of all remaining alignments is analyzed.
 *          For primary alignments, also the split read information is analyzed.
 */
void detect_junctions_in_alignment_file(const std::filesystem::path & alignment_file_path,
                                        const std::filesystem::path & insertion_file_path,
                                        const std::vector<uint8_t> methods,
                                        const clustering_methods clustering_method);
