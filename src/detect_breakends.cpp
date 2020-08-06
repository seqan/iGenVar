#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers

#include "detect_breakends/junction_detection.hpp"

struct cmd_arguments
{
    std::filesystem::path alignment_file_path{};
    std::filesystem::path insertion_file_path{};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "David Heller";
    parser.info.short_description = "Detect junctions in a read alignment file";
    parser.info.version = "0.0.1";
    parser.add_positional_option(args.alignment_file_path, "Input read alignments in SAM or BAM format.",
                                 seqan3::input_file_validator{{"sam", "bam"}} );
    parser.add_positional_option(args.insertion_file_path, "Output file for insertion alleles",
                                 seqan3::output_file_validator{{"fa", "fasta"}} );
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"detectJunctions", argc, argv};    // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);
    try
    {
        myparser.parse();                                               // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                   // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n";       // customise your error message
        return -1;
    }
    detect_junctions_in_alignment_file(args.alignment_file_path, args.insertion_file_path);

    return 0;
}
