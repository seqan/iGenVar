#pragma once

#include <forward_list>

#include <seqan3/io/sam_file/input.hpp>

/*! \brief Returns a list of pairs, where the first integer is a random chromosome and the second is a random position.
 *        The list is sorted by default.
 *  \param sample_size - The number of positions to generate.
 *  \param header - The header from the SAM file input.
 * 
 *  \return Returns a sorted forward_list of pairs of reference chromosome IDs and positions.
 */
std::forward_list<std::pair<int32_t, int32_t>> get_random_positions(uint64_t sample_size, seqan3::sam_file_header<std::deque<std::string>> & header)
{
    std::forward_list<std::pair<int32_t, int32_t>> output{};
    if (sample_size == 0)
    {
        output.push_front(std::pair<int32_t, int32_t>{-1, -1});
    }
    int32_t rand_chr{};
    int32_t rand_pos{};

    for (int i = 0; i < sample_size; ++i)
    {
        rand_chr = std::rand() % header.ref_ids().size();
        rand_pos = std::rand() % std::get<0>(header.ref_id_info[rand_chr]);

        output.push_front(std::pair<int32_t, int32_t>{rand_chr, rand_pos});
    }

    output.sort();
    return output;
}

/*! \brief Get the read length of a read from its cigar string.
 *  \param cigar - The cigar according to the seqan3 representation (a vector of cigar alphabet).
 *
 *  \return Returns the length of the cigar string.
 */

uint32_t get_read_length(std::vector<seqan3::cigar> & cigar)
{
    using seqan3::get;

    uint32_t result{0};
    for (auto e : cigar)
    {
        result += get<0>(e);
    }

    return result;
}
