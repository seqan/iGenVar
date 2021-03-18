#pragma once

#include "structures/aligned_segment.hpp"   // for struct AlignedSegment

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
 * \details The SA tag describes the alignments of a chimeric read and is like a small SAM within a SAM file:
 *          "SA:Z:(rname,pos,strand,CIGAR,mapQ,NM;)+"
 *          Each element (in parentheses) represents one alignment segment of the chimeric alignment formatted as
 *          a colon-delimited list.
 *          We add all segments to our candidate list `aligned_segments` and examine them in the following function
 *          `analyze_aligned_segments()`.
 *          For more information about this tag, see the
 *          ([Map Optional Fields Specification](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf)).
 */
void retrieve_aligned_segments(std::string const & sa_string, std::vector<AlignedSegment> & aligned_segments);

/*! \brief Build junctions out of aligned_segments.
 *
 * \param[in]       aligned_segments    - vector of [aligned_segments](\ref AlignedSegment)
 * \param[in, out]  junctions           - vector for storing junctions
 * \param[in, out]  query_sequence      - SEQ field of the SAM/BAM file
 * \param[in]       read_name           - QNAME field of the SAM/BAM file
 */
void analyze_aligned_segments(const std::vector<AlignedSegment> & aligned_segments,
                              std::vector<Junction> & junctions,
                              const seqan3::dna5_vector & query_sequence,
                              const std::string & read_name);

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
 * \param[in, out]  junctions   - vector for storing junctions
 */
void analyze_sa_tag(const std::string & query_name,
                    const seqan3::sam_flag & flag,
                    const std::string & ref_name,
                    const int32_t pos,
                    const uint8_t mapq,
                    std::vector<seqan3::cigar> const & cigar,
                    const seqan3::dna5_vector & seq,
                    const std::string & sa_tag,
                    std::vector<Junction> & junctions);
