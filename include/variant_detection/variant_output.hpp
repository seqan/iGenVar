#pragma once

#include <bio/var_io/writer.hpp>

#include "iGenVar.hpp"              // for cmd_arguments
#include "structures/cluster.hpp"   // for class Cluster

/*! \brief Gets the current time and transforms it in a nice readable way for the vcf header line filedate.
 *
 * \returns a time string in the format: YYYY-MM-DD HH:MM:SS.
 */
std::string transTime();

/*! \brief Writes a VCF header.
 *
 * \param[in] references_lengths - reference sequence dictionary parsed from \@SQ header lines
 * \param[in] sample_name        - Name of the sample for the vcf header line
 * \param[in, out] hdr           - reference sequence dictionary parsed from \@SQ header lines
 */
void write_header(std::map<std::string, int32_t> & references_lengths,
                  std::string sample_name,
                  bio::var_io::header & hdr);

/*! \brief Detects genomic variants from junction clusters and saves them in a record vector.
 *
 * \param[in] cluster       - input junction cluster
 * \param[in] args          - command line arguments:\n
 *                            **args.min_var_length** - minimum length of variants to detect
 *                                                      (expected to be non-negative) - *default: 30 bp* \n
 *                            **args.max_var_length** - maximum length of variants to detect
 *                                                      (expected to be non-negative) - *default: 1,000,000 bp*\n
 *                            **args.max_tol_inserted_length** - longest tolerated inserted sequence at non-INS SV types
 *                                                               (expected to be non-negative) - *default: 5 bp*\n
 *                            **args.min_qual** - minimum quality (amount of supporting reads) of a structural variant
 *                                                (expected to be non-negative) - *default: 1 supporting read*\n
 * \param[in, out] found_SV - will set to true, if an SV was found
 * \param[in, out] record   - vector of SV records
 *
 * \details Extracts genomic variants from given junction clusters.
 *          The quality of an SV is estimated based on the size of the cluster
 *          (i.e. the number of reads supporting the SV).
 */
void write_record(Cluster const & cluster,
                  cmd_arguments const & args,
                  bool & found_SV,
                  bio::var_io::default_record<> & record);

/*! \brief Detects genomic variants from junction clusters and prints them in output file in VCF format.
 *
 * \param[in] references_lengths - reference sequence dictionary parsed from \@SQ header lines
 * \param[in] clusters           - input junction clusters
 * \param[in] args               - command line arguments:\n
 *                                 **args.vcf_sample_name - Name of the sample for the vcf header line*\n
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
void find_and_output_variants(std::map<std::string, int32_t> & references_lengths,
                              std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::filesystem::path const & output_file_path);
