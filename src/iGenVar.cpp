#include <seqan3/argument_parser/argument_parser.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/range/views/get.hpp>

#include "variant_detection/variant_detection.hpp"
#include "variant_detection/validator.hpp"            // for class EnumValidator

// Specialise a mapping from an identifying string to the respective value of your type methods. With the help of this
// function, you're able to call ./detect_breackends with -m 1 and -m cigar_string and get the same result.
auto enumeration_names(detection_methods)
{
    return std::unordered_map<std::string, detection_methods>{{"0", detection_methods::cigar_string},
                                                              {"cigar_string", detection_methods::cigar_string},
                                                              {"1", detection_methods::split_read},
                                                              {"split_read", detection_methods::split_read},
                                                              {"2", detection_methods::read_pairs},
                                                              {"read_pairs", detection_methods::read_pairs},
                                                              {"3", detection_methods::read_depth},
                                                              {"read_depth", detection_methods::read_depth}};
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
    std::filesystem::path output_file_path{};
    std::vector<detection_methods> methods{cigar_string, split_read, read_pairs, read_depth};   // default: all methods
    clustering_methods clustering_method{simple_clustering};                                    // default: simple clustering method
    refinement_methods refinement_method{no_refinement};                                        // default: no refinement
    uint64_t min_var_length = 30;
    uint64_t sample_size{100000};                                                               // default is to sample 100000
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Lydia Buntrock, David Heller, Joshua Kim";
    parser.info.app_name = "iGenVar";
    parser.info.man_page_title = "Short and Long Read SV Caller";
    parser.info.short_description = "Detect genomic variants in a read alignment file";
    parser.info.version = "0.0.1";
    parser.info.date = "04-03-2021";    // last update
    parser.info.email = "lydia.buntrock@fu-berlin.de";
    parser.info.long_copyright = "long_copyright";
    parser.info.short_copyright = "short_copyright";
    parser.info.url = "https://github.com/seqan/iGenVar/";

    // Validatiors:
    seqan3::value_list_validator detection_method_validator {
        (seqan3::enumeration_names<detection_methods> | std::views::values)
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
    parser.add_option(args.output_file_path, 'o', "output",
                      "The path of the vcf output file. If no path is given, will output to standard output.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create,
                                                               {"vcf"}});

    // Options - Methods:
    parser.add_option(args.methods, 'm', "method", "Choose the detection method(s) to be used.",
                      seqan3::option_spec::advanced, detection_method_validator);
    parser.add_option(args.clustering_method, 'c', "clustering_method", "Choose the clustering method to be used.",
                      seqan3::option_spec::advanced, clustering_method_validator);
    parser.add_option(args.refinement_method, 'r', "refinement_method", "Choose the refinement method to be used.",
                      seqan3::option_spec::advanced, refinement_method_validator);

    // Options - SV specifications:
    parser.add_option(args.min_var_length, 'l', "min_var_length",
                      "Specify what should be the minimum length of your SVs to be detected (default 30 bp).",
                      seqan3::option_spec::advanced);

    // Options - Miscellaneous
    parser.add_option(args.sample_size, 's', "sample_size",
                      "This value is used in sampling positions for average read depth, insert size, and read length. Setting to 0 will skip sampling.",
                      seqan3::option_spec::advanced);
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"iGenVar", argc, argv};    // initialise myparser
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

    // Check that method selection contains no duplicates.
    std::vector<detection_methods> unique_methods{args.methods};
    std::ranges::sort(unique_methods);
    unique_methods.erase(std::unique(unique_methods.begin(), unique_methods.end()), unique_methods.end());
    if (args.methods.size() > unique_methods.size())
    {
        seqan3::debug_stream << "[Error] The same detection method was selected multiple times.\n";
        seqan3::debug_stream << "Methods to be used: " << args.methods << '\n';
        return -1;
    }

    detect_variants_in_alignment_file(args.alignment_file_path,
                                      args.insertion_file_path,
                                      args.methods,
                                      args.clustering_method,
                                      args.refinement_method,
                                      args.min_var_length,
                                      args.output_file_path,
                                      args.sample_size);

    return 0;
}
