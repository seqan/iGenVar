#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include "iGenVar.hpp"              // for global variable gVerbose
#include "structures/breakend.hpp"  // for class Breakend
#include "structures/junction.hpp"  // for class Junction

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;

std::span<seqan3::dna5 const> detect_tandem_duplication([[maybe_unused]] seqan3::dna5_vector const & query_sequence,
                                                        [[maybe_unused]] int32_t length,
                                                        [[maybe_unused]] int32_t pos_ref,
                                                        [[maybe_unused]] int32_t pos_read,
                                                        std::span<seqan3::dna5 const> & inserted_bases,
                                                        [[maybe_unused]] int32_t const min_length,
                                                        [[maybe_unused]] int32_t & pos_start_dup_seq,
                                                        [[maybe_unused]] int32_t & pos_end_dup_seq,
                                                        [[maybe_unused]] size_t & tandem_dup_count)
{
    // Suffix Case:
    // auto [ suffix_sequence_match_score, length_of_single_dupl_sequence_1 ] = align_suffix_or_prefix(...);

    // Prefix Case:
    // auto [ prefix_sequence_match_score, length_of_single_dupl_sequence_2 ] = align_suffix_or_prefix(...);

    // int16_t length_of_single_dupl_sequence = std::min(length_of_single_dupl_sequence_1,
    //                                                   length_of_single_dupl_sequence_2);
    // tandem_dup_count = suffix_sequence_match_score / length_of_single_dupl_sequence     // duplicated part on ref
    //                  + inserted_bases.size() / length_of_single_dupl_sequence           // inserted duplications
    //                  + prefix_sequence_match_score / length_of_single_dupl_sequence;    // duplicated part on ref
    // If its a duplication instead of an insertion, we save the (possible multiple times) duplicated part as duplicated_bases
    std::span<seqan3::dna5 const> duplicated_bases{};
    if (tandem_dup_count > 0)
        duplicated_bases = (inserted_bases /*| seqan3::views::slice(0, length_of_single_dupl_sequence)*/);
    return duplicated_bases;
}

void analyze_cigar(std::string const & read_name,
                   std::string const & chromosome,
                   int32_t const query_start_pos,
                   std::vector<seqan3::cigar> & cigar_string,
                   seqan3::dna5_vector const & query_sequence,
                   std::vector<Junction> & junctions,
                   int32_t const min_length)
{
    // Step through CIGAR string and store current position in reference and read
    int32_t pos_ref = query_start_pos;
    int32_t pos_read = 0;

    for (seqan3::cigar & pair : cigar_string)
    {
        using seqan3::get;
        int32_t length = get<0>(pair);
        seqan3::cigar::operation operation = get<1>(pair);
        if (operation == 'M'_cigar_operation || operation == '='_cigar_operation || operation == 'X'_cigar_operation)
        {
            pos_ref += length;
            pos_read += length;
        }
        // -------------------------------------
        else if (operation == 'I'_cigar_operation) // I: Insertion (gap in the reference sequence)
        {
            if (length >= min_length)
            {
                // Insertions cause one junction from the insertion location to the next base
                std::span<seqan3::dna5 const> inserted_bases = query_sequence | seqan3::views::slice(pos_read,
                                                                                                     pos_read + length);
        // ---------------- DUP ----------------
                // Case - Duplication: The inserted bases include one or more copies of a duplicated sequence, with an
                //                     origin somewhere else. -> global alignment
                // ##ALT=<ID=DUP,Description="Duplication">
                // TODO (23.7.21, irallia): This is Part of the Epic #144. Do we need a reference sequence for this?
        // ---------------- DUP:TANDEM ----------------
                size_t tandem_dup_count = 0;
                int32_t pos_start_dup_seq{};
                int32_t pos_end_dup_seq{};
                std::span<seqan3::dna5 const> duplicated_bases = detect_tandem_duplication(query_sequence,
                                                                                           length,
                                                                                           pos_ref,
                                                                                           pos_read,
                                                                                           inserted_bases,
                                                                                           min_length,
                                                                                           pos_start_dup_seq,
                                                                                           pos_end_dup_seq,
                                                                                           tandem_dup_count);
                if (tandem_dup_count != 0)
                {
                    Junction new_junction{Breakend{chromosome, pos_start_dup_seq, strand::forward},
                                          Breakend{chromosome, pos_end_dup_seq, strand::forward},
                                          duplicated_bases,
                                          tandem_dup_count,
                                          read_name};
                    if (gVerbose)
                        seqan3::debug_stream << "DUP:TANDEM: " << new_junction << "\n";
                    junctions.push_back(std::move(new_junction));
                } else {
        // ---------------- INS ----------------
                    Junction new_junction{Breakend{chromosome, pos_ref - 1, strand::forward},
                                          Breakend{chromosome, pos_ref, strand::forward},
                                          inserted_bases,
                                          0,
                                          read_name};
                    if (gVerbose)
                        seqan3::debug_stream << "INS: " << new_junction << "\n";
                    junctions.push_back(std::move(new_junction));
                }
            }
            pos_read += length;
        }
        // ---------------- DEL ----------------
        else if (operation == 'D'_cigar_operation) // D: Deletion (gap in the read)
        {
            if (length >= min_length)
            {
                // Deletions cause one junction from its start to its end
                Junction new_junction{Breakend{chromosome, pos_ref - 1, strand::forward},
                                      Breakend{chromosome, pos_ref + length, strand::forward},
                                      ""_dna5,
                                      0,
                                      read_name};
                if (gVerbose)
                    seqan3::debug_stream << "DEL: " << new_junction << "\n";
                junctions.push_back(std::move(new_junction));
            }
            pos_ref += length;
        }
        else if (operation == 'S'_cigar_operation)
        {
            pos_read += length;
        }
        else // other possible cigar operations: H, N, P
        {
            // seqan3::debug_stream << "Unhandled operation " << operation << std::endl;
        }

    }
}
