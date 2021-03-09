#pragma once

#include <ostream>

#include <seqan3/std/filesystem>

#include "structures/cluster.hpp"     // for class Cluster


/*! \brief Detects genomic variants from junction clusters and prints them to output stream in VCF format.
 *
 * \param clusters input junction clusters
 * \param out_stream output stream
 *
 * \details Extracts genomic variants from given junction clusters.
 */
void find_and_output_variant(std::vector<Cluster> const & clusters, std::ostream & out_stream);


/*! \brief Detects genomic variants from junction clusters and prints them in output file in VCF format.
 *
 * \param clusters input junction clusters
 * \param output_file_path output file path
 *
 * \details Extracts genomic variants from given junction clusters.
 */
//!\overload
void find_and_output_variant(std::vector<Cluster> const & clusters,
                             std::filesystem::path const & output_file_path);
