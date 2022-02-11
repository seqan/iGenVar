#pragma once

#include "structures/junction.hpp"  // for class Junction

/*! \brief This function checks if the inserted bases are tandem duplicated.
 *
 * \param[in] config                         - configuration for a semi-gobal alignment
 * \param[in] min_length                     - minimum length of variants to detect (default 30 bp,
 * \param[in] sequence                       - suffix or prefix sequence
 * \param[in] inserted_bases                 - the inserted bases of the possible duplication
 * \param[in] is_suffix                      - true: suffix case, false: prefix case
 * param[out] match_score                    - if we found a matching duplication, this value represents the maximal
 *                                             amount of matches with the reference and thus the length of the existing
 *                                             duplicated part on the reference. If there is no duplication, this value
 *                                             is 0.
 * param[out] length_of_single_dupl_sequence - greatest common divisor of length of inserted sequence and length of
 *                                             maching part -> length of a single duplicated sequence
 * \returns std::tie(match_score, length_of_single_dupl_sequence) - a tuple of the resulting values
 *
 * \details In this function the inserted bases are recursively aligned segment by segment until it has been proven that
 *          it is a real duplication.
 *          For simplicity, we only consider the suffix case here, since the prefix case works the same way:
 *
 *          suffix_sequence AAAACCGCGTAGCGGGGCGGG
 *                                     ||||||||||
 *          inserted_bases             GCGGGGCGGGGCGGG  -> unmatched_inserted_bases: GCGGG
 *                                                      -> match_score = 10
 *                                                      -> length_of_single_dupl_sequence = gcd(15, 10) = 5
 *
 *          suffix_sequence AAAACCGCGTAGCGGGGCGGG
 *                                          |||||
 *          unmatched_inserted_bases        GCGGG
 *
 *          -> match_score = 5, length_of_single_dupl_sequence = gcd(15, 5) = 5
 */
std::tuple<size_t, size_t> align_suffix_or_prefix(auto const & config,
                                                  int32_t const min_length,
                                                  std::span<const seqan3::dna5> & sequence,
                                                  std::span<const seqan3::dna5> & inserted_bases,
                                                  bool is_suffix);

/*! \brief This function checks if the inserted bases are tandem duplicated.
 *
 * \param[in]       query_sequence      - SEQ field of the SAM/BAM file
 * \param[in]       length              - length of inserted part, given by the CIGAR
 * \param[in]       pos_ref             - position of the inserted part in the ref (current position)
 * \param[in]       pos_read            - position of the inserted part in the read (current position)
 * \param[in]       inserted_bases      - the inserted bases of the possible duplication
 * \param[in]       min_length          - minimum length of variants to detect (default 30 bp,
 *                                                                              expected to be non-negative)
 * \param[in, out]  pos_start_dup_seq   - start position of the duplicated seq (excluding itself)
 * \param[in, out]  pos_end_dup_seq     - end position of the duplicated seq (including itself)
 * \param[in, out]  tandem_dup_count    - the number of tandem copies of the inserted sequence
 * \returns         duplicated_bases    - the duplicated bases of the duplication
 *
 * \details If the inserted bases include one or more copies of a duplicated sequence, which is suffix or prefix of
 *          another copy, we have a Tandem Duplication. To check this, we have to do a semi-global alignment.
 *
 *          Case 1: The duplication (insertion) comes after the matched sequence. Thus we need to check if the inserted
 *                  sequence matches (partly) the suffix sequence. The length of the matching part yields to the amount
 *                  of copies, thus we can calculate the tandem_dup_count and the inserted_bases.
 *
 *                  ref AAAACCGCGTAGCGGG----------TACGTAACGGTACG
 *                        ||||||||||||||          |||||||| -> inserted sequence: GCGGGGCGGG
 *                  read  AACCGCGTAGCGGGGCGGGGCGGGTACGTAAC
 *
 *                  suffix_sequence AAAACCGCGTAGCGGG       -> free_end_gaps_sequence1_leading{true},
 *                                             |||||          free_end_gaps_sequence1_trailing{false}
 *                  inserted_bases             GCGGGGCGGG  -> free_end_gaps_sequence2_leading{false},
 *                                                            free_end_gaps_sequence2_trailing{true}
 *                  -> tandem_dup_count = 3, duplicated_bases = GCGGG
 *
 *          Case 2: The duplication (insertion) comes before the matched sequence.
 *                  ref AAAACCGCGTA----------GCGGGTACGTAACGGTACG
 *                        |||||||||          |||||||||||||  -> inserted sequence: GCGGGGCGGG
 *                  read  AACCGCGTAGCGGGGCGGGGCGGGTACGTAAC
 *
 *                  prefix_sequence     GCGGGTACGTAACGGTACG -> free_end_gaps_sequence1_leading{false},
 *                                      |||||                  free_end_gaps_sequence1_trailing{true}
 *                  inserted_bases GCGGGGCGGG               -> free_end_gaps_sequence2_leading{true},
 *                                                             free_end_gaps_sequence2_trailing{false}
 *                  -> tandem_dup_count = 3, duplicated_bases = GCGGG
 * \see For some complex examples see detection_test.cpp.
 */
std::span<seqan3::dna5 const> detect_tandem_duplication(seqan3::dna5_vector const & query_sequence,
                                                        int32_t length,
                                                        int32_t pos_ref,
                                                        int32_t pos_read,
                                                        std::span<seqan3::dna5 const> & inserted_bases,
                                                        int32_t const min_length,
                                                        int32_t & pos_start_dup_seq,
                                                        int32_t & pos_end_dup_seq,
                                                        size_t & tandem_dup_count);

/*! \brief This function steps through the CIGAR string and stores junctions with their position in reference and read.
 *
 * \param[in]       read_name       - QNAME field of the SAM/BAM file
 * \param[in]       chromosome      - RNAME field of the SAM/BAM file
 * \param[in]       query_start_pos - POS field of the SAM/BAM file
 * \param[in]       cigar_string    - CIGAR field of the SAM/BAM file
 * \param[in]       query_sequence  - SEQ field of the SAM/BAM file
 * \param[in, out]  junctions       - vector for storing junctions
 * \param[in]       min_length      - minimum length of variants to detect (default 30 bp, expected to be non-negative)
 *
 * \details This function steps through the CIGAR string and stores junctions with their position in reference and read.
 *          We distinguish 4 cases of CIGAR operation characters:
 *          1. M, =, X: For matches and mismatches (Alignment column containing two letters. This could contain two
 *                      different letters (mismatch) or two identical letters), we step through ref and read.
 *          2. I:       For insertions (gap in the reference sequence) -> Insertions cause two junctions ( (1) from the
 *                      reference to the read and (2) back from the read to the reference ).
 *          3. D:       For deletions (gap in the query sequence) -> Deletions cause one junction from its start to its
 *                      end.
 *          4. S:       For soft clipped letters we step through the read. These are segments of the query sequence that
 *                      do not appear in the alignment. The full-length query sequence is given in the SEQ field of the
 *                      SAM record.
 *
 *          Other CIGAR operations: H, N, P are skipped (H: hard clipping sequences are not present in the SEQ, N:
 *          skipped region representing an intron, P: padding consumes neither the query nor the reference).
 *          The junctions found are stored in the given `junctions` vector.
 *          For more information about CIGAR operations see the
 *          [Map Format Specification](https://samtools.github.io/hts-specs/SAMv1.pdf#page=8) (last access 09.04.2021).
 */
void analyze_cigar(std::string const & read_name,
                   std::string const & chromosome,
                   int32_t const query_start_pos,
                   std::vector<seqan3::cigar> & cigar_string,
                   seqan3::dna5_vector const & query_sequence,
                   std::vector<Junction> & junctions,
                   int32_t const min_length);
