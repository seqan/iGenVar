#pragma once

/*! \brief This function steps through the CIGAR string and stores junctions with their position in reference and read.
 *
 * \param read_name         QNAME field of the SAM/BAM file
 * \param chromosome        RNAME field of the SAM/BAM file
 * \param query_start_pos   POS field of the SAM/BAM file
 * \param cigar_string      CIGAR field of the SAM/BAM file
 * \param query_sequence    SEQ field of the SAM/BAM file
 * \param junctions         vector for storing junctions
 * \param min_length        minimum length of variants to detect (default 30 bp)
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
 *          Other CIGAR operations: H, N, P are skipped (H: hard clipping sequences are not present in the SEQ, N:
 *          skipped region representing an intron, P: padding consumes neither the query nor the reference).
 *          The junctions found are stored in the given `junctions` vector.
 *          For more information see the
 *          ([Map Format Specification](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf)) page 8.
 */
void analyze_cigar(const std::string & read_name,
                   const std::string chromosome,
                   const int32_t query_start_pos,
                   std::vector<seqan3::cigar> & cigar_string,
                   const seqan3::dna5_vector & query_sequence,
                   std::vector<Junction> & junctions,
                   uint64_t const min_length);
