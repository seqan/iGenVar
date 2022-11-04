#pragma once

#include <filesystem>
#include <map>
#include <set>
#include <string>
#include <vector>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "structures/junction.hpp"

//! \brief Stores information from the genome file.
struct Genome
{
    std::vector<seqan3::dna5_vector> seqs; //!< The sequences.
    std::vector<std::string> names; //!< The names of the sequences.
};

/*!
 * \brief Read multiple genome sequences.
 * \param genome_filename The filename where to read from.
 * \return a genome object with sequences and names.
 */
Genome read_genome(std::filesystem::path const & genome_filename);

//! \brief Stores information about read activity, i.e. diversions from a genome.
struct Activity
{
    std::vector<std::vector<unsigned>> values; //!< A measure for read activity.
    std::vector<int32_t> refmap; //!< A mapping from the genome sequence id to the reference id.
};

/*!
 * \brief Read a sam/bam file and retrieve vectors of activity for each genome position.
 * \param[in] reads_filename The file with reads (sam/bam).
 * \param[in] genome A genome object with sequences and names.
 * \param[in] min_var_length The length above which an indel/SNP is considered a variant.
 * \return an activity object with a value vector for each genome sequence.
 */
Activity analyze_activity(std::filesystem::path const & reads_filename, Genome const & genome, uint64_t min_var_length);

//! \brief The type for storing active regions (an interval on the genome where activity is high).
typedef std::vector<std::pair<int, int>> ActiveRegions;

/*!
 * \brief Extract regions of high activity from an activity profile.
 * \param[in] activity - The activity values for the reference genomes.
 * \param[in] window_width - The size of the sliding window (high values lead to larger regions).
 * \param[in] act_threshold - The threshold when a window is considered active.
 * \return lists of region coordinates: the first is the start position, the second is one-after-the-end position.
 */
std::vector<ActiveRegions> get_active_regions(Activity activity, size_t window_width = 25, size_t act_threshold = 2);

/*!
 * \brief Store a SNP or indel that has been found in a haplotype sequence.
 * \param[out] junctions The set of junctions where the SNP is added.
 * \param[in] bufR The deleted substring in the reference sequence.
 * \param[in] bufH The inserted substring in the haplotype sequence (alternative).
 * \param[in] pos The genome position where the SNP ends.
 * \param[in] ref_name The name of the corresponding reference sequence.
 * \param[in] qual A quality value for the reliability of this change.
 *
 * \details
 * In the case of an insertion the reference substring is empty,
 * in the case of a deletion the haplotype substring is empty.
 * If none of them is empty, we have found a substitution.
 */
void store_snp(std::set<Junction> & junctions,
               seqan3::dna5_vector & bufR,
               seqan3::dna5_vector & bufH,
               int pos,
               std::string const & ref_name,
               float qual);

/*!
 * \brief Detect Single Nucleotide Polymorphisms (SNPs) and short deletions and insertions.
 * \param[out] junctions - The resulting SNPs are appended to this set.
 * \param[out] references_lengths - A dictionary for storing information about reference sequences.
 * \param[in] genome_filename - The file path where to find the genome reference sequence.
 * \param[in] reads_filename - The file path where to find the sequenced reads.
 * \param[in] min_var_length - The length above which an indel/SNP is considered a variant.
 */
void detect_snp_and_indel(std::set<Junction> & junctions,
                          std::map<std::string, int32_t> & references_lengths,
                          std::filesystem::path const & genome_filename,
                          std::filesystem::path const & reads_filename,
                          uint64_t min_var_length);
