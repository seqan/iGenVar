#pragma once

#include <seqan3/std/filesystem>

/*!
 * \brief Detect Single Nucleotide Polymorphisms (SNPs) and short deletions and insertions.
 * \param[in] reads_filename - The file path where to find the sequenced reads.
 */
void detect_snp_and_indel(std::filesystem::path const & reads_filename);
