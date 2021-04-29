#pragma once

#include <ostream>

#include <seqan3/std/filesystem>

#include "iGenVar.hpp"              // for cmd_arguments
#include "structures/cluster.hpp"   // for class Cluster


/*! \brief Detects genomic variants from junction clusters and prints them to output stream in VCF format.
 *
 * \param[in]       clusters    - input junction clusters
 * \param[in]       args        - command line arguments:\n
 *                                **args.min_var_length** - minimum length of variants to detect - *default: 30 bp*
 * \param[in, out]  out_stream  - output stream
 *
 * \details Extracts genomic variants from given junction clusters.
 */
void find_and_output_variants(std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::ostream & out_stream);


/*! \brief Detects genomic variants from junction clusters and prints them in output file in VCF format.
 *
 * \param[in] clusters          - input junction clusters
 * \param[in] args              - command line arguments:\n
 *                                **args.min_var_length** - minimum length of variants to detect - *default: 30 bp*
 * \param[in] output_file_path  - output file path
 *
 * \details Extracts genomic variants from given junction clusters.
 */
//!\overload
void find_and_output_variants(std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::filesystem::path const & output_file_path);
