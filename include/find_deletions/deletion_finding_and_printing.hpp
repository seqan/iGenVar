#pragma once

#include <seqan3/std/filesystem>
#include <fstream>
#include <iostream>
#include <vector>

#include "junction.hpp"

void find_and_print_deletions(const std::filesystem::path & junction_file_path);
