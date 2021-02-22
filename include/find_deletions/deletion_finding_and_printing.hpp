#pragma once

#include <seqan3/std/filesystem>
#include <ostream>

/*! \brief Detects deletions out of the junction file.
 *
 * \param junction_file_path input junction file
 * \param out_stream output stream
 *
 * \details Extracts deletions out of given breakends / junctions.
 */
void find_and_print_deletions(std::filesystem::path const & junction_file_path, std::ostream & out_stream);

//!\overload
void find_and_print_deletions(std::filesystem::path const & junction_file_path,
                              std::filesystem::path const & output_file_path);
