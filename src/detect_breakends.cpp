#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/get.hpp>

#include "detect_breakends/junction_detection.hpp"
#include "detect_breakends/validator.hpp"            // for class EnumValidator

// Specialise a mapping from an identifying string to the respective value of your type methods. With the help of this
// function, you're able to call ./detect_breackends with -m 1 and -m cigar_string and get the same result.
auto enumeration_names(detecting_methods)
{
    return std::unordered_map<std::string, detecting_methods>{{"0", detecting_methods::cigar_string},
                                                              {"cigar_string", detecting_methods::cigar_string},
                                                              {"1", detecting_methods::split_read},
                                                              {"split_read", detecting_methods::split_read},
                                                              {"2", detecting_methods::read_pairs},
                                                              {"read_pairs", detecting_methods::read_pairs},
                                                              {"3", detecting_methods::read_depth},
                                                              {"read_depth", detecting_methods::read_depth}};
};

// Specialise a mapping from an identifying string to the respective value of your type clustering_methods. With the
// help of this function, you're able to call ./detect_breackends with -c 0 and -c simple_clustering and get the same
// result.
auto enumeration_names(clustering_methods)
{
    return std::unordered_map<std::string,
                              clustering_methods>{{"0", clustering_methods::simple_clustering},
                                                  {"simple_clustering",
                                                   clustering_methods::simple_clustering},
                                                  {"1", clustering_methods::hierarchical_clustering},
                                                  {"hierarchical_clustering",
                                                   clustering_methods::hierarchical_clustering},
                                                  {"2", clustering_methods::self_balancing_binary_tree},
                                                  {"self_balancing_binary_tree",
                                                   clustering_methods::self_balancing_binary_tree},
                                                  {"3", clustering_methods::candidate_selection_based_on_voting},
                                                  {"candidate_selection_based_on_voting",
                                                   clustering_methods::candidate_selection_based_on_voting}};
};

// Specialise a mapping from an identifying string to the respective value of your type refinement_methods. With the
// help of this function, you're able to call ./detect_breackends with -r 0 and -r no_refinement and get the same
// result.
auto enumeration_names(refinement_methods)
{
    return std::unordered_map<std::string, refinement_methods>{{"0", refinement_methods::no_refinement},
                                                               {"no_refinement",
                                                                refinement_methods::no_refinement},
                                                               {"1", refinement_methods::sViper_refinement_method},
                                                               {"sViper_refinement_method",
                                                                refinement_methods::sViper_refinement_method},
                                                               {"2", refinement_methods::sVirl_refinement_method},
                                                               {"sVirl_refinement_method",
                                                                refinement_methods::sVirl_refinement_method}};
};

struct cmd_arguments
{
    std::filesystem::path alignment_file_path{};
    std::filesystem::path insertion_file_path{};
    std::vector<detecting_methods> methods{cigar_string, split_read, read_pairs, read_depth};   // default: all methods
    clustering_methods clustering_method{simple_clustering};    // default: simple clustering method
    refinement_methods refinement_method{no_refinement};        // default: no refinement
    uint64_t min_var_length = 30;
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "David Heller & Lydia Buntrock";
    parser.info.app_name = "iGenVar";
    parser.info.man_page_title = "Short and long Read SV Caller";
    parser.info.short_description = "Detect junctions in a read alignment file";
    parser.info.version = "0.0.1";
    parser.info.date = "19-01-2021";    // last update
    parser.info.email = "lydia.buntrock@fu-berlin.de";
    parser.info.long_copyright = "long_copyright";
    parser.info.short_copyright = "short_copyright";
    parser.info.url = "https://github.com/seqan/iGenVar/";

    // Validatiors:
    seqan3::value_list_validator detecting_method_validator {
        (seqan3::enumeration_names<detecting_methods> | std::views::values)
    };
    // ToDo (Lydia): Should get solved with solving https://github.com/seqan/iGenVar/issues/78
    // EnumValidator<detecting_methods> detecting_method_validator{seqan3::enumeration_names<detecting_methods>
    //                                                             | std::views::values};
    EnumValidator<clustering_methods> clustering_method_validator{seqan3::enumeration_names<clustering_methods>
                                                                  | std::views::values};
    EnumValidator<refinement_methods> refinement_method_validator{seqan3::enumeration_names<refinement_methods>
                                                                  | std::views::values};

    // Options - Input / Output:
    parser.add_positional_option(args.alignment_file_path, "Input read alignments in SAM or BAM format.",
                                 seqan3::input_file_validator{{"sam", "bam"}} );
    parser.add_positional_option(args.insertion_file_path, "Output file for insertion alleles.",
                                 seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create,
                                                               {"fa", "fasta"}} );

    // Options - Methods:
    parser.add_option(args.methods, 'm', "method", "Choose the detecting method(s) to be used.",
                      seqan3::option_spec::advanced, detecting_method_validator);
    parser.add_option(args.clustering_method, 'c', "clustering_method", "Choose the clustering method to be used.",
                      seqan3::option_spec::advanced, clustering_method_validator);
    parser.add_option(args.refinement_method, 'r', "refinement_method", "Choose the refinement method to be used.",
                      seqan3::option_spec::advanced, refinement_method_validator);

    // Options - SV specifications:
    parser.add_option(args.min_var_length, 'l', "min_var_length",
                      "Specify what should be the minimum length of your SVs to be detected (default 30 bp).",
                      seqan3::option_spec::advanced);
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"detectJunctions", argc, argv};    // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);

    // Parse the given arguments and catch possible errors.
    try
    {
        myparser.parse();                                               // trigger command line parsing
    }
    catch (seqan3::argument_parser_error const & ext)                   // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << '\n';       // customise your error message
        return -1;
    }

    seqan3::debug_stream << "Methods to be used: " << args.methods << '\n';

    detect_junctions_in_alignment_file(args.alignment_file_path,
                                       args.insertion_file_path,
                                       args.methods,
                                       args.clustering_method,
                                       args.refinement_method,
                                       args.min_var_length);

    return 0;
}
