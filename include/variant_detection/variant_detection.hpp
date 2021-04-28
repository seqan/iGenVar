#pragma once

#include <seqan3/std/filesystem>    // for filesystem
#include <vector>

#include "iGenVar.hpp"              // for struct cmd_arguments
#include "structures/junction.hpp"  // for class Junction

/*! \brief Detects junctions between distant genomic positions by analyzing a short read alignment file (sam/bam). The
 *         detected junctions are stored in a vector.
 *
 * \param[in, out]  junctions - a vector of junctions
 * \param[in]       args - command line arguments:\n
 *                         **args.alignment_short_reads_file_path** - short reads input file, path to the sam/bam file\n
 *                         **args.methods** - list of methods for detecting junctions
 *                            (0: cigar_string, 1: split_read, 2: read_pairs, 3: read_depth) - *default: all methods*\n
 * 
 *
 * \details Detects junctions from the CIGAR strings and supplementary alignment tags of read alignment records.
 *          We filter unmapped alignments, secondary alignments, duplicates and alignments with low mapping quality.
 *          Then, the CIGAR string of all remaining alignments is analyzed.
 *          For primary alignments, also the split read information is analyzed.
 */
void detect_junctions_in_short_reads_sam_file(std::vector<Junction> & junctions,
                                              cmd_arguments const & args);

/*! \brief Detects junctions between distant genomic positions by analyzing a long read alignment file (sam/bam). The
 *         detected junctions are stored in a vector.
 *
 * \param[in, out]  junctions - a vector of junctions
 * \param[in]       args - command line arguments:\n
 *                         **args.alignment_long_reads_file_path** - long reads input file, path to the sam/bam file\n
 *                         **args.methods** - list of methods for detecting junctions
 *                            (0: cigar_string, 1: split_read, 2: read_pairs, 3: read_depth) - *default: all methods*\n
 *                         **args.min_var_length** - minimum length of variants to detect - *default: 30 bp*\n
 *                         **args.max_overlap** - maximum overlap between alignment segments - *default: 10 bp*
 *
 *
 * \details Detects junctions from the CIGAR strings and supplementary alignment tags of read alignment records.
 *          We filter unmapped alignments, secondary alignments, duplicates and alignments with low mapping quality.
 *          Then, the CIGAR string of all remaining alignments is analyzed.
 *          For primary alignments, also the split read information is analyzed.
 */
void detect_junctions_in_long_reads_sam_file(std::vector<Junction> & junctions,
                                             cmd_arguments const & args);
