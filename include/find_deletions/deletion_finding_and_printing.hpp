#pragma once

#include <seqan3/std/filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "junction.hpp"

/*! \brief Detects deletions out of the junction file.
 *
 * \param junction_file_path input junction file
 *
 * \details Extracts deletions out of given breakends / junctions.
 */
void find_and_print_deletions(const std::filesystem::path & junction_file_path);
