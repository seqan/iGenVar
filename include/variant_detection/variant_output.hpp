#pragma once

#include <filesystem>
#include <ostream>

#include "iGenVar.hpp"              // for cmd_arguments
#include "structures/cluster.hpp"   // for class Cluster


/*! \brief Detects genomic variants from junction clusters and prints them to output stream in VCF format.
 *
 * \param[in] references_lengths - reference sequence dictionary parsed from \@SQ header lines
 * \param[in] clusters           - input junction clusters
 * \param[in] args               - command line arguments:\n
 *                                 **args.min_var_length** - minimum length of variants to detect
 *                                                           (expected to be non-negative) - *default: 30 bp* \n
 *                                 **args.max_var_length** - maximum length of variants to detect
 *                                                           (expected to be non-negative) - *default: 1,000,000 bp*\n
 *                                 **args.max_tol_inserted_length** - longest tolerated inserted sequence at non-INS SV
 *                                                                    types (expected to be non-negative)
 *                                                                  - *default: 5 bp*\n
 *                                 **args.min_qual** - minimum quality (amount of supporting reads) of a structural
 *                                                     variant (expected to be non-negative)
 *                                                   - *default: 1 supporting read*\n
 * \param[in, out] out_stream    - output stream
 *
 * \details Extracts genomic variants from given junction clusters.
 *          The quality of an SV is estimated based on the size of the cluster
 *          (i.e. the number of reads supporting the SV).
 */
void find_and_output_variants(std::map<std::string, int32_t> & references_lengths,
                              std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::ostream & out_stream);


/*! \brief Detects genomic variants from junction clusters and prints them in output file in VCF format.
 *
 * \param[in] references_lengths - reference sequence dictionary parsed from \@SQ header lines
 * \param[in] clusters           - input junction clusters
 * \param[in] args               - command line arguments:\n
 *                                 **args.min_var_length** - minimum length of variants to detect
 *                                                           (expected to be non-negative) - *default: 30 bp* \n
 *                                 **args.max_var_length** - maximum length of variants to detect
 *                                                           (expected to be non-negative) - *default: 1,000,000 bp*\n
 *                                 **args.max_tol_inserted_length** - longest tolerated inserted sequence at non-INS SV
 *                                                                    types (expected to be non-negative)
 *                                                                  - *default: 5 bp*\n
 *                                 **args.min_qual** - minimum quality (amount of supporting reads) of a structural
 *                                                     variant (expected to be non-negative)
 *                                                   - *default: 1 supporting read*\n
 ** \param[in] output_file_path  - output file path
 *
 * \details Extracts genomic variants from given junction clusters.
 *          The quality of an SV is estimated based on the size of the cluster
 *          (i.e. the number of reads supporting the SV).
 */
//!\overload
void find_and_output_variants(std::map<std::string, int32_t> & references_lengths,
                              std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::filesystem::path const & output_file_path);
