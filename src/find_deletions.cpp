#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers

#include "find_deletions/deletion_finding_and_printing.hpp"

struct cmd_arguments
{
    std::filesystem::path junction_file_path{};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "David Heller";
    parser.info.short_description = "Find deletions";
    parser.info.version = "0.0.1";
    parser.add_positional_option(args.junction_file_path, "Input junctions tab-separated format.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"partitionJunctions", argc, argv}; // initialise myparser
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
    find_and_print_deletions(args.junction_file_path);


    return 0;
}
