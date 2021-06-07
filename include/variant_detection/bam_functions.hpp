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
