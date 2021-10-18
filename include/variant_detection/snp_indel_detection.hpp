#pragma once

#include <seqan3/std/filesystem>

/*!
 * \brief Detect Single Nucleotide Polymorphisms (SNPs) and short deletions and insertions.
 * \param[in] reads_filename - The file path where to find the sequenced reads.
 * \param[in] min_var_length - The length above which an indel/SNP is considered a variant.
 */
void detect_snp_and_indel(std::filesystem::path const & reads_filename, uint64_t min_var_length);
