#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include <seqan3/core/debug_stream.hpp>     // our custom output stream
#include <seqan3/std/filesystem>            // use std::filesystem::path
#include <seqan3/io/sequence_file/all.hpp>  // FASTA support
#include <seqan3/io/alignment_file/all.hpp> // SAM/BAM support

using namespace seqan3;

void run_program(std::filesystem::path & alignment_file_path, std::filesystem::path & reference_file_path)
{
    alignment_file_input alignment_file{alignment_file_path, fields<field::ID, field::SEQ, field::FLAG>{}};
    sequence_file_input reference_file{reference_file_path};

    for (auto & [id, seq, flag] : alignment_file)
    {
        debug_stream << id << seq << flag << std::endl;
    }
}

struct cmd_arguments
{
    std::filesystem::path alignment_file_path{};
    std::filesystem::path reference_file_path{};
};

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "David Heller";
    parser.info.short_description = "Detect breakends in a read alignment file";
    parser.info.version = "0.0.1";
    parser.add_positional_option(args.alignment_file_path, "Read alignments in SAM or BAM format.",
                                 input_file_validator{{"sam", "bam"}} );
    parser.add_positional_option(args.reference_file_path, "Reference genome in FASTA format",
                                 input_file_validator{{"fa", "fasta"}} );
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
    run_program(args.alignment_file_path, args.reference_file_path);
    return 0;
}
