#pragma once

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/exception.hpp>
#include <seqan3/io/sam_file/output.hpp>
#include <seqan3/utility/char_operations/predicate.hpp>


/*!\details
 * From the SAM/BAM Format Specification
 * | Bit  | Bit   | Description                                                        | Note                                               |
 * |:----:|:-----:|:------------------------------------------------------------------:|:--------------------------------------------------:|
 * |    1 | 0x1   | template having multiple segments in sequencing                    | necessary for read pair method                     |
 * |    2 | 0x2   | each segment properly aligned according to the aligner             | aligned as expected, probably not containing SVs *)|
 * |    4 | 0x4   | segment unmapped                                                   | will be skipped                                    |
 * |    8 | 0x8   | next segment in the template unmapped                              | mate unmapped, maybe an insertion? *)              |
 * |   16 | 0x10  | SEQ being reverse complemented                                     | important for split reads and paired-reads         |
 * |   32 | 0x20  | SEQ of the next segment in the template being reverse complemented | *)                                                 |
 * |   64 | 0x40  | the first segment in the template                                  | *)                                                 |
 * |  128 | 0x80  | the last segment in the template                                   | *)                                                 |
 * |  256 | 0x100 | secondary alignment                                                | will be skipped                                    |
 * |  512 | 0x200 | not passing filters, such as platform/vendor quality controls      | *)                                                 |
 * | 1024 | 0x400 | PCR or optical duplicate                                           | will be skipped                                    |
 * | 2048 | 0x800 | supplementary alignment                                            | will be skipped (primary alignments only)          |
 * *) not checked yet
 * \see [Map Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf#page=7)
 * (last access 09.04.2021).
 *
 */
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

inline constexpr bool hasFlagMultiple(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_MULTIPLE) == BAM_FLAG_MULTIPLE;
}

// inline constexpr bool hasFlagProperMappedPairedReads(seqan3::sam_flag const & flag)
// {
//     return (static_cast<uint16_t>(flag) & BAM_FLAG_ALL_PROPER) == BAM_FLAG_ALL_PROPER;
// }

inline constexpr bool hasFlagUnmapped(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_UNMAPPED) == BAM_FLAG_UNMAPPED;
}

/*!\brief This flag is important for split reads and paired reads. For the latter, 0x10 and 0x20 are compared. Usually
 *        you would expect them to be different and if they are not it could be because of an inversion.
 */
//TODO (irallia): Implement the use of this flag in the paired read method.
inline constexpr bool hasFlagReverseComplement(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_RC) == BAM_FLAG_RC;
}

inline constexpr bool hasFlagSecondary(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_SECONDARY) == BAM_FLAG_SECONDARY;
}

inline constexpr bool hasFlagSupplementary(seqan3::sam_flag const & flag)
{
    return (static_cast<uint16_t>(flag) & BAM_FLAG_SUPPLEMENTARY) == BAM_FLAG_SUPPLEMENTARY;
}

inline constexpr bool hasFlagDuplicate(seqan3::sam_flag const & flag)
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
 * \tparam cigar_input_type - the type of a single pass input view over the cigar string; must model
 *                            std::ranges::input_range
 * \param[in, out] cigar_input - the single pass input view over the cigar string to parse
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
inline constexpr std::tuple<std::vector<seqan3::cigar>, int32_t, int32_t> parse_cigar(cigar_input_type && cigar_input)
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
        auto buff_end = (std::ranges::copy(cigar_view | std::views::take_while(seqan3::is_digit), buffer.data())).out;
        cigar_operation = *std::ranges::begin(cigar_view);
        std::ranges::next(std::ranges::begin(cigar_view));

        if (std::from_chars(buffer.begin(), buff_end, cigar_count).ec != std::errc{})
            throw seqan3::format_error{"Corrupted cigar string encountered"};

        update_alignment_lengths(ref_length, seq_length, cigar_operation, cigar_count);
        operations.emplace_back(cigar_count, seqan3::cigar::operation{}.assign_char(cigar_operation));
    }

    return {operations, ref_length, seq_length};
}
