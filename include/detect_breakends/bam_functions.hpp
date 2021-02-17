#pragma once

#include <seqan3/range/views/all.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>
#include <seqan3/io/exception.hpp>

enum BamFlags
{
    BAM_FLAG_MULTIPLE      = 0x0001,
    BAM_FLAG_ALL_PROPER    = 0x0002,
    BAM_FLAG_UNMAPPED      = 0x0004,
    BAM_FLAG_NEXT_UNMAPPED = 0x0008,
    BAM_FLAG_RC            = 0x0010,
    BAM_FLAG_NEXT_RC       = 0x0020,
    BAM_FLAG_FIRST         = 0x0040,
    BAM_FLAG_LAST          = 0x0080,
    BAM_FLAG_SECONDARY     = 0x0100,
    BAM_FLAG_QC_NO_PASS    = 0x0200,
    BAM_FLAG_DUPLICATE     = 0x0400,
    BAM_FLAG_SUPPLEMENTARY = 0x0800
};

inline bool
hasFlagUnmapped(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_UNMAPPED) == BAM_FLAG_UNMAPPED;
}

inline bool
hasFlagReverseComplement(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_RC) == BAM_FLAG_RC;
}

inline bool
hasFlagSecondary(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_SECONDARY) == BAM_FLAG_SECONDARY;
}

inline bool
hasFlagSupplementary(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_SUPPLEMENTARY) == BAM_FLAG_SUPPLEMENTARY;
}

inline bool
hasFlagDuplicate(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_DUPLICATE) == BAM_FLAG_DUPLICATE;
}


 /*! \brief Updates the sequence lengths by `cigar_count` depending on the cigar operation `op`.
  * \param[in, out]  ref_length  The reference sequence's length.
  * \param[in, out]  seq_length  The query sequence's length.
  * \param[in]       op          The cigar operation.
  * \param[in]       cigar_count The cigar count value to add to the length depending on the cigar operation.
  */
static void update_alignment_lengths(int32_t & ref_length, int32_t & seq_length, char op, uint32_t cigar_count)
  {
switch (op)
{
    case 'M': case '=': case 'X': ref_length += cigar_count, seq_length += cigar_count; break;
    case 'D': case 'N':           ref_length += cigar_count; break;
    case 'I': case 'S':           seq_length += cigar_count; break;
    case 'H': case 'P':           break; // no op (hard-clipping or padding does not increase either length)
    default:                      throw seqan3::format_error{"Illegal cigar operation: " + std::string{op}};
}
};


/*! \brief Parses a cigar string into a vector of operation-count pairs (e.g. (M, 3)).
 * \tparam cigar_input_type The type of a single pass input view over the cigar string; must model
 *                          std::ranges::input_range.
 * \param[in]  cigar_input  The single pass input view over the cigar string to parse.
 *
 * \returns A tuple of size three containing (1) std::vector over seqan3::cigar, that describes
 *          the alignment, (2) the aligned reference length, (3) the aligned query sequence length.
 *
 * \details
 *
 * For example, the view over the cigar string "1H4M1D2M2S" will return
 * `{[(H,1), (M,4), (D,1), (M,2), (S,2)], 7, 6}`.
 */
template <typename cigar_input_type>
inline std::tuple<std::vector<seqan3::cigar>, int32_t, int32_t> parse_cigar(cigar_input_type && cigar_input)
{
    std::vector<seqan3::cigar> operations{};
    std::array<char, 20> buffer{}; // buffer to parse numbers with from_chars. Biggest number should fit in uint64_t
    char cigar_operation{};
    uint32_t cigar_count{};
    int32_t ref_length{}, seq_length{}; // length of aligned part for ref and query

    // transform input into a single input view if it isn't already
    auto cigar_view = cigar_input | seqan3::views::single_pass_input;

    // parse the rest of the cigar
    // -------------------------------------------------------------------------------------------------------------
    while (std::ranges::begin(cigar_view) != std::ranges::end(cigar_view)) // until stream is not empty
    {
        auto buff_end = (std::ranges::copy(cigar_view | seqan3::views::take_until_or_throw(!seqan3::is_digit), buffer.data())).out;
        cigar_operation = *std::ranges::begin(cigar_view);
        std::ranges::next(std::ranges::begin(cigar_view));

        if (std::from_chars(buffer.begin(), buff_end, cigar_count).ec != std::errc{})
            throw seqan3::format_error{"Corrupted cigar string encountered"};

        update_alignment_lengths(ref_length, seq_length, cigar_operation, cigar_count);
        operations.emplace_back(cigar_count, seqan3::cigar_op{}.assign_char(cigar_operation));
    }

    return {operations, ref_length, seq_length};
}
