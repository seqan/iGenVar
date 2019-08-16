#include <string>
#include <vector>

#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>     // our custom output stream
#include <seqan3/std/filesystem>            // use std::filesystem::path
#include <seqan3/io/sequence_file/all.hpp>  // FASTA support
#include <seqan3/io/alignment_file/all.hpp> // SAM/BAM support
#include <seqan3/alphabet/nucleotide/dna5.hpp>

#include "bam_functions.h"

using namespace seqan3;

class Breakend {
    bool in_reference1;
    int32_t chromosome1;
    int32_t position1;
    bool forward1;
    bool in_reference2;
    int32_t chromosome2;
    int32_t position2;
    bool forward2;
    public:
        Breakend (bool, int32_t, int32_t, bool, bool, int32_t, int32_t, bool);
        void print_vcf_entry ();
} ;


Breakend::Breakend (bool in_ref1, int32_t chr1, int32_t pos1, bool fwd1, bool in_ref2, int32_t chr2, int32_t pos2, bool fwd2)
{
    in_reference1 = in_ref1;
    chromosome1 = chr1;
    position1 = pos1;
    forward1 = fwd1;
    in_reference2 = in_ref2;
    chromosome2 = chr2;
    position2 = pos2;
    forward2 = fwd2;
}

void Breakend::print_vcf_entry()
{
    printf("%d\t%d\t%s\t%d\t%d\t%s\n", chromosome1, position1, forward1 ? "-->" : "<--", chromosome2, position2, forward2 ? "-->" : "<--");
}


void analyze_cigar(std::vector<cigar> & cigar_string, std::vector<Breakend> & breakends, std::vector<dna5_vector> & insertions, int32_t chromosome, int32_t query_start_pos, dna5_vector & query_sequence, int32_t min_length, sequence_file_output<> & insertion_file)
{
    // Step through CIGAR string and store current position in reference and read
    int32_t pos_ref = query_start_pos;
    int32_t pos_read = 0;

    // Stores the index of the current read in the insertion allele output file (or -1 if current read has not been added yet)
    int32_t insertion_allele_id {-1};

    for (cigar & pair : cigar_string)
    {
        int32_t length = get<0>(pair);
        cigar_op operation = get<1>(pair);
        if (operation == 'M'_cigar_op || operation == '='_cigar_op || operation == 'X'_cigar_op)
        {
            pos_ref += length;
            pos_read += length;
        }
        else if (operation == 'I'_cigar_op)
        {
            if (length > min_length)
            {
                // debug_stream << "Found insertion of length " << length << " at position " << pos_ref << '\n';
                if (insertion_allele_id < 0)
                {
                    insertion_allele_id = insertions.size();
                    std::string insertion_allele_name{"allele_" + std::to_string(insertion_allele_id)};
                    insertion_file.emplace_back(query_sequence, insertion_allele_name);
                    insertions.push_back(query_sequence);
                }
                // Insertions cause two breakends ( (1) from the reference to the read and (2) back from the read to the reference )
                Breakend new_breakend1 = Breakend {true, chromosome, query_start_pos, true, false, insertion_allele_id, pos_read, true};
                new_breakend1.print_vcf_entry();
                breakends.push_back(new_breakend1);
                Breakend new_breakend2 = Breakend {false, insertion_allele_id, pos_read + length, true, true, chromosome, query_start_pos, true};
                new_breakend2.print_vcf_entry();
                breakends.push_back(new_breakend2);
            }
            pos_read += length;
        }
        else if (operation == 'D'_cigar_op)
        {
            if (length > min_length)
            {
                // debug_stream << "Found deletion of length " << length << " at position " << pos_ref << '\n';
                // Deletions cause one breakend from its start to its end
                Breakend new_breakend = Breakend {true, chromosome, query_start_pos, true, true, chromosome, query_start_pos + length, true};
                new_breakend.print_vcf_entry();
                breakends.push_back(new_breakend);
            }
            pos_ref += length;
        }
        else if (operation == 'S'_cigar_op)
        {
            pos_read += length;
        }       
    }
}

void detect_breakends_in_alignment_file(std::filesystem::path & alignment_file_path, std::filesystem::path & insertion_file_path)
{
    // Open input alignment file
    using my_fields = fields<field::REF_ID, field::REF_OFFSET, field::FLAG, field::MAPQ, field::CIGAR, field::SEQ>;
    alignment_file_input alignment_file{alignment_file_path, my_fields{}};

    // Open output file for insertion alleles
    sequence_file_output insertion_file{insertion_file_path};

    // Store breakends, insertion_alleles and number of good alignments
    std::vector<Breakend> breakends {};
    std::vector<dna5_vector> insertion_alleles {};
    uint16_t num_good = 0;

    for (auto & rec : alignment_file)
    {
        int32_t chr = get<field::REF_ID>(rec).value_or(0);
        int32_t pos = get<field::REF_OFFSET>(rec).value_or(0);
        auto flag = get<field::FLAG>(rec);
        auto mapq = get<field::MAPQ>(rec);
        auto cigar = get<field::CIGAR>(rec);
        auto seq = get<field::SEQ>(rec);

        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20)
        {
            // debug_stream << "Skipped flag " << flag << std::endl;
        }
        else
        {
            assert(pos >= 0);
            analyze_cigar(cigar, breakends, insertion_alleles, chr, pos, seq, 30, insertion_file);
            num_good++;
            if (num_good % 1000 == 0)
            {
                debug_stream << num_good << " good alignments" << std::endl;
            }
        }
    }

    debug_stream << "Done. Found " << breakends.size() << " breakends." << '\n';
}

struct cmd_arguments
{
    std::filesystem::path alignment_file_path{};
    std::filesystem::path insertion_file_path{};
};

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "David Heller";
    parser.info.short_description = "Detect breakends in a read alignment file";
    parser.info.version = "0.0.1";
    parser.add_positional_option(args.alignment_file_path, "Input read alignments in SAM or BAM format.",
                                 input_file_validator{{"sam", "bam"}} );
    parser.add_positional_option(args.insertion_file_path, "Output file for insertion alleles",
                                 output_file_validator{{"fa", "fasta"}} );
}

int main(int argc, char ** argv)
{
    argument_parser myparser{"detectBreakends", argc, argv};        // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);
    try
    {
        myparser.parse();                                          // trigger command line parsing
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        debug_stream << "[Error] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
    detect_breakends_in_alignment_file(args.alignment_file_path, args.insertion_file_path);
    return 0;
}