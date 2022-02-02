#pragma once

#include <filesystem>

#include <seqan3/alphabet/cigar/cigar.hpp>

/*!
 * \brief Detect Single Nucleotide Polymorphisms (SNPs) and short deletions and insertions.
 * \param[in] reads_filename - The file path where to find the sequenced reads.
 * \param[in] min_var_length - The length above which an indel/SNP is considered a variant.
 */
void detect_snp_and_indel(std::filesystem::path const & reads_filename, uint64_t min_var_length);

/*!
 * \brief Extract activity from SAM records by counting indels and soft clips.
 * \param[in,out] activity       - The activity values for one reference genome.
 * \param[in]     cigar_sequence - The cigar string of the SAM record.
 * \param[in]     min_var_length - The length above which an indel/SNP is considered a variant.
 * \param[in]     ref_pos        - The start position of the alignment in the genome.
 */
void update_activity_for_record(std::vector<unsigned> & activity,
                                std::vector<seqan3::cigar> const & cigar_sequence,
                                uint64_t min_var_length,
                                int32_t ref_pos);

/*!
 * \brief Extract active regions from activity profile.
 * \param[in] activity - The activity values for one reference genome.
 * \return a list of position intervals, where the activity is high.
 */
std::vector<std::pair<size_t, size_t>> active_regions(std::vector<unsigned> const & activity);
