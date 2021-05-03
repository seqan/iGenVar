#pragma once

#include <seqan3/std/filesystem>    // for filesystem
#include <vector>

#include "method_enums.hpp"         // for enum detection_methods, clustering_methods and refinement_methods
#include "structures/junction.hpp"  // for class Junction

/*! \brief Detects junctions between distant genomic positions by analyzing a short read alignment file (sam/bam). The
 *         detected junctions are stored in a vector.
 *
 * \param[in, out]  junctions - a vector of junctions
 * \param[in]       alignment_short_reads_file_path - short reads input file, path to the sam/bam file
 * \param[in]       methods - list of methods for detecting junctions (0: cigar_string,
 *                                                                     1: split_read,
 *                                                                     2: read_pairs,
 *                                                                     3: read_depth)
 * \param[in]       min_var_length - minimum length of variants to detect (default 30 bp)
 *
 * \details Detects junctions from the CIGAR strings and supplementary alignment tags of read alignment records.
 *          We filter unmapped alignments, secondary alignments, duplicates and alignments with low mapping quality.
 *          Then, the CIGAR string of all remaining alignments is analyzed.
 *          For primary alignments, also the split read information is analyzed.
 */
void detect_junctions_in_short_reads_sam_file(std::vector<Junction> & junctions,
                                              std::filesystem::path const & alignment_short_reads_file_path,
                                              std::vector<detection_methods> const & methods,
                                              uint64_t const min_var_length);

/*! \brief Detects junctions between distant genomic positions by analyzing a long read alignment file (sam/bam). The
 *         detected junctions are stored in a vector.
 *
 * \param[in, out]  junctions - a vector of junctions
 * \param[in]       alignment_long_reads_file_path - long reads input file, path to the sam/bam file
 * \param[in]       methods - list of methods for detecting junctions (0: cigar_string,
 *                                                                     1: split_read,
 *                                                                     2: read_pairs,
 *                                                                     3: read_depth)
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
                                             uint64_t const min_var_length);
