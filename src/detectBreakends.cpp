#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>     // our custom output stream
#include <seqan3/std/filesystem>            // use std::filesystem::path
#include <seqan3/io/sequence_file/all.hpp>  // FASTA support
#include <seqan3/io/alignment_file/all.hpp> // SAM/BAM support
#include "bam_functions.h"

#include <string>
#include <vector>

using namespace seqan3;

class Breakend {
    int32_t chromosome1;
    int32_t position1;
    bool forward1;
    int32_t chromosome2;
    int32_t position2;
    bool forward2;
    public:
        Breakend (int32_t, int32_t, bool, int32_t, int32_t, bool);
        void print_vcf_entry ();
} ;


Breakend::Breakend (int32_t chr1, int32_t pos1, bool fwd1, int32_t chr2, int32_t pos2, bool fwd2)
{
    chromosome1 = chr1;
    position1 = pos1;
    forward1 = fwd1;
    chromosome2 = chr2;
    position2 = pos2;
    forward2 = fwd2;
}

void Breakend::print_vcf_entry()
{
    printf("%d\t%d\t%s\t%d\t%d\t%s\n", chromosome1, position1, forward1 ? "-->" : "<--", chromosome2, position2, forward2 ? "-->" : "<--");
}


void analyze_cigar(std::vector<cigar> & cigar_string, std::vector<Breakend> & breakends, int32_t chromosome, int32_t query_start_pos, int32_t min_length)
{
    int32_t pos_ref = 0;
    int32_t pos_read = 0;
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
                debug_stream << "Found insertion of length " << length << " at position " << pos_ref << '\n';
            }
            pos_read += length;
        }
        else if (operation == 'D'_cigar_op)
        {
            if (length > min_length)
            {
                debug_stream << "Found deletion of length " << length << " at position " << pos_ref << '\n';
                Breakend new_breakend {chromosome, query_start_pos + pos_ref, true, chromosome, query_start_pos + pos_ref + length, true};
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


void run_program(std::filesystem::path & alignment_file_path)
{
    using my_fields = fields<field::REF_ID, field::REF_OFFSET, field::FLAG, field::MAPQ, field::CIGAR>;
    alignment_file_input alignment_file{alignment_file_path, my_fields{}};

    uint16_t num_good = 0;
    std::vector<Breakend> breakends {};
    for (auto & rec : alignment_file)
    {
        int32_t chr = get<field::REF_ID>(rec).value_or(0);
        int32_t pos = get<field::REF_OFFSET>(rec).value_or(0);
        auto flag = get<field::FLAG>(rec);
        auto mapq = get<field::MAPQ>(rec);
        auto cigar = get<field::CIGAR>(rec);

        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20)
        {
            // debug_stream << "Skipped flag " << flag << std::endl;
        }
        else
        {
            assert(pos >= 0);
            analyze_cigar(cigar, breakends, chr, pos, 30);
            num_good++;
            if (num_good % 1000 == 0)
            {
                std::cout << num_good << " good alignments" << std::endl;
            }
        }
    }

    debug_stream << "Found " << breakends.size() << " breakends." << '\n';
    for (Breakend b : breakends)
    {
        b.print_vcf_entry();
    }
}

struct cmd_arguments
{
    std::filesystem::path alignment_file_path{};
    // std::filesystem::path reference_file_path{};
};

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "David Heller";
    parser.info.short_description = "Detect breakends in a read alignment file";
    parser.info.version = "0.0.1";
    parser.add_positional_option(args.alignment_file_path, "Read alignments in SAM or BAM format.",
                                 input_file_validator{{"sam", "bam"}} );
    // parser.add_positional_option(args.reference_file_path, "Reference genome in FASTA format",
    //                             input_file_validator{{"fa", "fasta"}} );
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
    run_program(args.alignment_file_path);
    return 0;
}