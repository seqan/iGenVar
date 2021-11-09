#pragma once

#include <seqan3/std/filesystem>    // for filesystem
#include <map>
#include <vector>

#include "iGenVar.hpp"              // for struct cmd_arguments
#include "structures/junction.hpp"  // for class Junction

#include "bamit/all.hpp"

/*! \brief Reads the header of the input file. Checks if input file is sorted and reads the reference sequence
 *         dictionary. Stores the reference sequence lengths in parameter `reference_lengths` and returns the list of
 *         reference sequences.
 *
 * \param[in, out]  alignment_file - short or long reads input file
 * \param[in, out]  references_lengths - reference sequence dictionary parsed from \@SQ header lines
 */
std::deque<std::string> read_header_information(auto & alignment_file,
                                                std::map<std::string, int32_t> & references_lengths);

/*! \brief Attempts to load a bamit index having the same name as a given input file, with ".bit" appended at the end.
 *         If this file does not exist, it will create the index itself and save it to that file.
 *
 * \param[in] input_path - the path to the short/long read alignment file.
 * \return A bamit::index tree.
 */
std::vector<std::unique_ptr<bamit::IntervalNode>> load_or_create_index(std::filesystem::path const & input_path);

/*! \brief Detects junctions between distant genomic positions by analyzing a short read alignment file (sam/bam). The
 *         detected junctions are stored in a vector.
 *
 * \param[in, out]  junctions - a vector of junctions
 * \param[in, out]  references_lengths - reference sequence dictionary parsed from \@SQ header lines
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
                                              std::map<std::string, int32_t> & references_lengths,
                                              cmd_arguments const & args);

/*! \brief Detects junctions between distant genomic positions by analyzing a long read alignment file (sam/bam). The
 *         detected junctions are stored in a vector.
 *
 * \param[in, out]  junctions - a vector of junctions
 * \param[in, out]  references_lengths - reference sequence dictionary parsed from \@SQ header lines
 * \param[in]       args - command line arguments:\n
 *                         **args.alignment_long_reads_file_path** - long reads input file, path to the sam/bam file\n
 *                         **args.methods** - list of methods for detecting junctions
 *                            (0: cigar_string, 1: split_read, 2: read_pairs, 3: read_depth) - *default: all methods*\n
 *                         **args.min_var_length** - minimum length of variants to detect
 *                            (expected to be non-negative) - *default: 30 bp*\n
 *                         **args.max_overlap** - maximum overlap between alignment segments
 *                            (expected to be non-negative) - *default: 10 bp*
 *
 *
 * \details Detects junctions from the CIGAR strings and supplementary alignment tags of read alignment records.
 *          We filter unmapped alignments, secondary alignments, duplicates and alignments with low mapping quality.
 *          Then, the CIGAR string of all remaining alignments is analyzed.
 *          For primary alignments, also the split read information is analyzed.
 */
void detect_junctions_in_long_reads_sam_file(std::vector<Junction> & junctions,
                                             std::map<std::string, int32_t> & references_lengths,
                                             cmd_arguments const & args);
