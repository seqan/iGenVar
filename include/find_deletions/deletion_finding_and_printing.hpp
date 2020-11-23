#pragma once

#include <seqan3/std/filesystem>
#include <fstream>
#include <vector>

#include "junction.hpp"

/*! \brief Detects deletions out of the junction file.
 *
 * \param junction_file_path input junction file
 * \param output_file_path output vcf file
 *
 * \details Extracts deletions out of given breakends / junctions.
 */
void find_and_print_deletions(std::filesystem::path const & junction_file_path, std::filesystem::path const & output_file_path);
