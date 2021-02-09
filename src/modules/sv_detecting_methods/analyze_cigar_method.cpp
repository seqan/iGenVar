#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/alphabet/cigar/cigar_op.hpp>
#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sequence_file/output.hpp>

#include "structures/junction.hpp"  // for class Junction
#include "structures/breakend.hpp"  // for class Breakend

using seqan3::operator""_cigar_op;

void analyze_cigar(const std::string & read_name,
                   const std::string chromosome,
                   const int32_t query_start_pos,
                   std::vector<seqan3::cigar> & cigar_string,
                   const seqan3::dna5_vector & query_sequence,
                   std::vector<Junction> & junctions,
                   std::vector<seqan3::dna5_vector> & insertions,
                   int32_t min_length,
                   seqan3::sequence_file_output<> & insertion_file)
{
    // Step through CIGAR string and store current position in reference and read
    int32_t pos_ref = query_start_pos;
    int32_t pos_read = 0;

    // Stores the index of the current read in the insertion allele output file (or -1 if current read has not been added yet)
    int32_t insertion_allele_id {-1};

    for (seqan3::cigar & pair : cigar_string)
    {
        using seqan3::get;
        int32_t length = get<0>(pair);
        seqan3::cigar_op operation = get<1>(pair);
        if (operation == 'M'_cigar_op || operation == '='_cigar_op || operation == 'X'_cigar_op)
        {
            pos_ref += length;
            pos_read += length;
        }
        else if (operation == 'I'_cigar_op) // I: Insertion (gap in the reference sequence)
        {
            if (length >= min_length)
            {
                if (insertion_allele_id < 0)
                {
                    insertion_allele_id = insertions.size();
                    std::string insertion_allele_name{"allele_" + std::to_string(insertion_allele_id)};
                    insertion_file.emplace_back(query_sequence, insertion_allele_name);
                    insertions.push_back(query_sequence);
                }
                // Insertions cause two junctions ( (1) from the reference to the read and (2) back from the read to the reference )
                Junction new_junction1{Breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                       Breakend{std::to_string(insertion_allele_id),
                                                pos_read,
                                                strand::forward,
                                                sequence_type::read},
                                       read_name};
                seqan3::debug_stream << "INS1: " << new_junction1 << "\n";
                junctions.push_back(std::move(new_junction1));
                Junction new_junction2{Breakend{std::to_string(insertion_allele_id),
                                                pos_read + length,
                                                strand::forward,
                                                sequence_type::read},
                                       Breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                       read_name};
                seqan3::debug_stream << "INS2: " << new_junction2 << "\n";
                junctions.push_back(std::move(new_junction2));
            }
            pos_read += length;
        }
        else if (operation == 'D'_cigar_op)
        {
            if (length >= min_length)
            {
                // Deletions cause one junction from its start to its end
                Junction new_junction{Breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                      Breakend{chromosome, pos_ref + length, strand::forward, sequence_type::reference},
                                      read_name};
                seqan3::debug_stream << "DEL: " << new_junction << "\n";
                junctions.push_back(std::move(new_junction));
            }
            pos_ref += length;
        }
        else if (operation == 'S'_cigar_op)
        {
            pos_read += length;
        }
        else // other possible cigar operations: H, N, P
        {
            // seqan3::debug_stream << "Unhandled operation " << operation << std::endl;
        }

    }
}
