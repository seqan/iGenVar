#pragma once

#include "iGenVar.hpp"                          // for struct cmd_arguments
#include "structures/aligned_segment.hpp"       // for struct AlignedSegment
#include "structures/junction.hpp"              // for class Junction
#include "variant_detection/bam_functions.hpp"  // for seqan3::sam_flag and hasFlag* functions

/*! \brief Splits a string by a given delimiter and stores substrings in a given container.
 *
 * \param[in]       str     - string to split
 * \param[in, out]  cont    - container for the splitted substrings
 * \param[in]       delim   - delimiter
 */
template <class Container>
void split_string(std::string const & str, Container & cont, char const delim = ' ');

/*! \brief Parse the SA tag from the SAM/BAM alignment of a chimeric/split-aligned read. Build
 *         [aligned_segments](\ref AlignedSegment), one for each alignment segment of the read.
 *
 * \param[in]       sa_string           - "SA" tag string
 * \param[in, out]  aligned_segments    - vector of [aligned_segments](\ref AlignedSegment)
 *
 * \details The SA tag describes the alignments of a chimeric read and is like a small SAM within a SAM
 *          file:
 *
 *          `"SA:Z:(rname,pos,strand,CIGAR,mapQ,NM;)+"`
 *
 *          Each element (in parentheses) represents one alignment segment of the chimeric alignment formatted as
 *          a colon-delimited list.
 *          We add all segments to our candidate list `aligned_segments` and examine them in the following function
 *          `analyze_aligned_segments()`.
 *
 *          For more information about this tag, see the
 *          [Map Optional Fields Specification](https://samtools.github.io/hts-specs/SAMtags.pdf)
 *          (last access 09.04.2021).
 */
void retrieve_aligned_segments(std::string const & sa_string, std::vector<AlignedSegment> & aligned_segments);

/*! \brief Calculates the end respectively start position of two consecutive mates (current and next
 *         [AlignedSegment](\ref AlignedSegment)) depending on their orientation and corrects the position of the second
 *         mate, if they are overlapping.
 *
 * \param[in] current          - current [AlignedSegment](\ref AlignedSegment)
 * \param[in] next             - consecutive next [AlignedSegment](\ref AlignedSegment)
 * \param[in] distance_on_read - distance between the consecutive aligned segments on the read
 *
 * \returns std::pair(mate1_pos, mate2_pos) - a pair of the resulting mate positions
 *
 * \details \verbatim
    Case 1: --current--> --next-->
    Case 2: --current--> <--next--
    Case 3: <--current-- --next-->
    Case 4: <--current-- <--next--
                       | |
               mate1_pos mate2_pos

    overlapping Case: --current-->         distance_on_read < 0
                                --next-->
                                 ||
                         mate1_pos mate2_pos
    \endverbatim
 */
std::pair<int32_t, int32_t> get_mate_positions(AlignedSegment current, AlignedSegment next, int32_t distance_on_read);

/*! \brief Build junctions out of aligned_segments.
 *
 * \param[in]       aligned_segments    - vector of [aligned_segments](\ref AlignedSegment)
 * \param[in, out]  junctions           - vector for storing junctions
 * \param[in, out]  query_sequence      - SEQ field of the SAM/BAM file
 * \param[in]       read_name           - QNAME field of the SAM/BAM file
 * \param[in]       min_length          - minimum length of variants to detect (expected to be non-negative)
 * \param[in]       max_overlap         - maximum overlap between alignment segments (expected to be non-negative)
 */
void analyze_aligned_segments(std::vector<AlignedSegment> const & aligned_segments,
                              std::vector<Junction> & junctions,
                              seqan3::dna5_vector const & query_sequence,
                              std::string const & read_name,
                              int32_t const min_length,
                              int32_t const max_overlap);

/*! \brief Parse the SA tag from the SAM/BAM alignment of a chimeric/split-aligned read. Build
 *         [aligned_segments](\ref AlignedSegment), one for each alignment segment of the read.
 *         Sort the alignment segments and analyze them to detect junctions.
 *
 * \param[in]       query_name  - QNAME field of the SAM/BAM file
 * \param[in]       flag        - FLAG field of the SAM/BAM file
 * \param[in]       ref_name    - RNAME field of the SAM/BAM file
 * \param[in]       pos         - POS field of the SAM/BAM file
 * \param[in]       mapq        - MAPQ field of the SAM/BAM file
 * \param[in]       cigar       - CIGAR field of the SAM/BAM file
 * \param[in]       seq         - SEQ field of the SAM/BAM file
 * \param[in]       sa_tag      - SA tag, one tag from the read of the SAM/BAM file
 * \param[in]       args        - command line arguments:\n
 *                                **args.min_var_length** - minimum length of variants to detect (expected to be non-negative)\n
 *                                **args.max_overlap** - maximum overlap between alignment segments (expected to be non-negative)
 * \param[in, out]  junctions   - vector for storing junctions
 */
void analyze_sa_tag(std::string const & query_name,
                    seqan3::sam_flag const & flag,
                    std::string const & ref_name,
                    int32_t const pos,
                    uint8_t const mapq,
                    std::vector<seqan3::cigar> const & cigar,
                    seqan3::dna5_vector const & seq,
                    std::string const & sa_tag,
                    cmd_arguments const & args,
                    std::vector<Junction> & junctions);
