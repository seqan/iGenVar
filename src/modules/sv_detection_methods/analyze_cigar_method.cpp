#include <limits>   // for std::numeric_limits
#include <numeric>  // for std::gcd (greatest common divisor)

// #include <seqan3/alignment/aligned_sequence/debug_stream_alignment.hpp> // for a nicer alignment view
#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_scoring_scheme.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alignment/scoring/nucleotide_scoring_scheme.hpp>
#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include "iGenVar.hpp"              // for global variable gVerbose
#include "structures/breakend.hpp"  // for class Breakend
#include "structures/junction.hpp"  // for class Junction

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;

std::tuple<size_t, size_t> align_suffix_or_prefix(auto const & config,
                                                  int32_t const min_length,
                                                  std::span<seqan3::dna5 const> & sequence,
                                                  std::span<seqan3::dna5 const> & inserted_bases,
                                                  bool is_suffix)
{
    // length of matching part of suffix or prefix sequence
    size_t match_score{};
    // The number of tandem copies of this junction.
    size_t length_of_single_dupl_sequence = std::numeric_limits<size_t>::max();

    auto results = seqan3::align_pairwise(std::tie(sequence, inserted_bases), config);
    auto & res = *results.begin();
    // TODO (irallia 17.8.21): The mismatches should give us the opportunity to allow a given amount of errors in the
    // duplication.
    size_t matches = res.score() % 100;
    size_t mismatches = (res.score() - matches) * (-1);
    // For the beginning we do not allow mistakes, we should change this later, see todo above.
    if (mismatches > 0)
        return std::tie(match_score, length_of_single_dupl_sequence);
    // found duplicated sequence in front of inserted sequence
    if (matches >= min_length)
    {
        // seqan3::debug_stream << "Resulting alignment:\n" << res.alignment() << '\n';
        match_score = matches;
        // The possible length of the single duplicated: greatest common divisor of length of inserted part and length
        // of maching part
        length_of_single_dupl_sequence = std::gcd(inserted_bases.size(), matches);
        if(matches != inserted_bases.size())
        {
            std::span<seqan3::dna5 const> unmatched_inserted_bases{};
            if (is_suffix)
                unmatched_inserted_bases = inserted_bases | seqan3::views::slice(match_score, inserted_bases.size());
            else
                unmatched_inserted_bases = inserted_bases | seqan3::views::slice(0, inserted_bases.size() - match_score);
            // The first length_of_single_dupl_sequence is already the greatest common divisor and therefore the next
            // recursively calculated one can be ignored.
            auto [ next_match_score,
                   next_length_of_single_dupl_sequence ] = align_suffix_or_prefix(config,
                                                                                  min_length,
                                                                                  sequence,
                                                                                  unmatched_inserted_bases,
                                                                                  is_suffix);
            // If the substring does not match, there is no real duplication.
            if (next_match_score == 0) {
                length_of_single_dupl_sequence = std::numeric_limits<size_t>::max();
                return std::tie(next_match_score, length_of_single_dupl_sequence);
            }
        }
    }
    // The first match_score calculated is automatically the maximum, so the recursively next one can be ignored.
    return std::tie(match_score, length_of_single_dupl_sequence);
}

std::span<seqan3::dna5 const> detect_tandem_duplication(seqan3::dna5_vector const & query_sequence,
                                                        int32_t length,
                                                        int32_t pos_ref,
                                                        int32_t pos_read,
                                                        std::span<seqan3::dna5 const> & inserted_bases,
                                                        int32_t const min_length,
                                                        int32_t & pos_start_dup_seq,
                                                        int32_t & pos_end_dup_seq,
                                                        size_t & tandem_dup_count)
{
    auto scoring_scheme = seqan3::align_cfg::scoring_scheme{
                              seqan3::nucleotide_scoring_scheme{seqan3::match_score{1},
                                                                seqan3::mismatch_score{-100}}};

    auto gap_scheme = seqan3::align_cfg::gap_cost_affine{seqan3::align_cfg::open_score{0},
                                                         seqan3::align_cfg::extension_score{-100}};
    // Suffix Case:
    auto suffix_config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                          seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                                          seqan3::align_cfg::free_end_gaps_sequence1_trailing{false},
                                                          seqan3::align_cfg::free_end_gaps_sequence2_trailing{true}} |
                         scoring_scheme | gap_scheme;

    std::span<seqan3::dna5 const> suffix_sequence = query_sequence | seqan3::views::slice(0, pos_read);
    auto [ suffix_sequence_match_score, length_of_single_dupl_sequence_1 ] = align_suffix_or_prefix(suffix_config,
                                                                                                    min_length,
                                                                                                    suffix_sequence,
                                                                                                    inserted_bases,
                                                                                                    true);
    // Prefix Case:
    auto prefix_config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{false},
                                                          seqan3::align_cfg::free_end_gaps_sequence2_leading{true},
                                                          seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                          seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}} |
                         scoring_scheme | gap_scheme;

    std::span<seqan3::dna5 const> prefix_sequence = query_sequence | seqan3::views::slice(pos_read + length,
                                                                                          query_sequence.size());
    auto [ prefix_sequence_match_score, length_of_single_dupl_sequence_2 ] = align_suffix_or_prefix(prefix_config,
                                                                                                    min_length,
                                                                                                    prefix_sequence,
                                                                                                    inserted_bases,
                                                                                                    false);

    pos_start_dup_seq = pos_ref - 1 - suffix_sequence_match_score;
    pos_end_dup_seq = pos_ref - 1 + prefix_sequence_match_score;

    int16_t length_of_single_dupl_sequence = std::min(length_of_single_dupl_sequence_1,
                                                      length_of_single_dupl_sequence_2);
    tandem_dup_count = suffix_sequence_match_score / length_of_single_dupl_sequence     // duplicated part on ref
                     + inserted_bases.size() / length_of_single_dupl_sequence           // inserted duplications
                     + prefix_sequence_match_score / length_of_single_dupl_sequence;    // duplicated part on ref
    // If its a duplication instead of an insertion, we save the (possible multiple times) duplicated part as duplicated_bases
    std::span<seqan3::dna5 const> duplicated_bases{};
    if (tandem_dup_count > 0)
        duplicated_bases = (inserted_bases | seqan3::views::slice(0, length_of_single_dupl_sequence));
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
                    {
                        seqan3::debug_stream << "DUP:TANDEM: " << new_junction << "\n";
                        seqan3::debug_stream << "\t\t\tduplicated sequence: " << duplicated_bases
                                             << " with " << tandem_dup_count << " duplications\n";
                    }
                    junctions.push_back(std::move(new_junction));
                } else {
        // ---------------- INS ----------------
                    Junction new_junction{Breakend{chromosome, pos_ref - 1, strand::forward},
                                          Breakend{chromosome, pos_ref, strand::forward},
                                          inserted_bases,
                                          0,
                                          read_name};
                    if (gVerbose)
                    {
                        seqan3::debug_stream << "INS: " << new_junction << "\n";
                        seqan3::debug_stream << "\t\t\tinserted sequence: " << inserted_bases << "\n";
                    }
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
