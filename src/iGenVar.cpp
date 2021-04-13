#include "iGenVar.hpp"

#include <seqan3/core/debug_stream.hpp>                     // for seqan3::debug_stream

#include "modules/clustering/simple_clustering_method.hpp"  // for the simple clustering method
#include "structures/cluster.hpp"                           // for class Cluster
#include "variant_detection/validator.hpp"                  // for class EnumValidator
#include "variant_detection/variant_detection.hpp"          // for detect_junctions_in_long_reads_sam_file()
#include "variant_detection/variant_output.hpp"             // for find_and_output_variants()

struct cmd_arguments
{
    std::filesystem::path alignment_short_reads_file_path{""};
    std::filesystem::path alignment_long_reads_file_path{""};
    std::filesystem::path output_file_path{};
    std::vector<detection_methods> methods{cigar_string, split_read, read_pairs, read_depth};   // default: all methods
    clustering_methods clustering_method{simple_clustering};                                    // default: simple clustering method
    refinement_methods refinement_method{no_refinement};                                        // default: no refinement
    uint64_t min_var_length = 30;
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Lydia Buntrock, David Heller, Joshua Kim";
    parser.info.app_name = "iGenVar";
    parser.info.man_page_title = "Short and Long Read SV Caller";
    parser.info.short_description = "Detect genomic variants in a read alignment file";
    parser.info.version = "0.0.2";
    parser.info.date = "30-03-2021";    // last update
    parser.info.email = "lydia.buntrock@fu-berlin.de";
    parser.info.long_copyright = "long_copyright";
    parser.info.short_copyright = "short_copyright";
    parser.info.url = "https://github.com/seqan/iGenVar/";

    // Validatiors:
    EnumValidator<detection_methods> detection_method_validator{seqan3::enumeration_names<detection_methods>
                                                                | std::views::values};
    EnumValidator<clustering_methods> clustering_method_validator{seqan3::enumeration_names<clustering_methods>
                                                                  | std::views::values};
    EnumValidator<refinement_methods> refinement_method_validator{seqan3::enumeration_names<refinement_methods>
                                                                  | std::views::values};

    // Options - Input / Output:
    parser.add_option(args.alignment_short_reads_file_path,
                      'i', "input_short_reads",
                      "Input short read alignments in SAM or BAM format (Illumina).",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"sam", "bam"}} );
    parser.add_option(args.alignment_long_reads_file_path,
                      'j', "input_long_reads",
                      "Input long read alignments in SAM or BAM format (PacBio, Oxford Nanopore, ...).",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator{{"sam", "bam"}} );
    parser.add_option(args.output_file_path, 'o', "output",
                      "The path of the vcf output file. If no path is given, will output to standard output.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"vcf"}});

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
}

void detect_variants_in_alignment_file(cmd_arguments const & args)
{
    // Store junctions
    std::vector<Junction> junctions{};

    // short reads
    // ToDo (Lydia): handle short reads
    if (args.alignment_short_reads_file_path != "")
    {
        seqan3::debug_stream << "Short reads are currently not supported.\n ";
        return;
    }

    // long reads
    if (args.alignment_long_reads_file_path == "")
    {
        seqan3::debug_stream << "No long reads were given (short reads are currently not supported).\n ";
        return;
    }
    detect_junctions_in_long_reads_sam_file(junctions,
                                            args.alignment_long_reads_file_path,
                                            args.methods,
                                            args.clustering_method,
                                            args.refinement_method,
                                            args.min_var_length);

    seqan3::debug_stream << "Start clustering...\n";

    std::vector<Cluster> clusters{};
    switch (args.clustering_method)
    {
        case 0: // simple_clustering
            simple_clustering_method(junctions, clusters);
            break;
        case 1: // hierarchical clustering
            seqan3::debug_stream << "The hierarchical clustering method is not yet implemented\n";
            break;
        case 2: // self-balancing_binary_tree,
            seqan3::debug_stream << "The self-balancing binary tree clustering method is not yet implemented\n";
            break;
        case 3: // candidate_selection_based_on_voting
            seqan3::debug_stream << "The candidate selection based on voting clustering method is not yet implemented\n";
            break;
    }

    seqan3::debug_stream << "Done with clustering. Found " << clusters.size() << " junction clusters.\n";

    switch (args.refinement_method)
    {
        case 0: // no refinement
            seqan3::debug_stream << "No refinement was selected.\n";
            break;
        case 1: // sViper_refinement_method
            seqan3::debug_stream << "The sViper refinement method is not yet implemented\n";
            break;
        case 2: // sVirl_refinement_method
            seqan3::debug_stream << "The sVirl refinement method is not yet implemented\n";
            break;
    }

    find_and_output_variants(clusters, args.output_file_path);
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

    // Check if we have at least one input file.
    if (args.alignment_short_reads_file_path == "" && args.alignment_long_reads_file_path == "")
    {
        seqan3::debug_stream << "[Error] You need to input at least one sam/bam file.\n"
                             << "Please use -i or -input_short_reads to pass a short read file "
                                "or -j or -input_long_reads for a long read file.\n";
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

    detect_variants_in_alignment_file(args);

    return 0;
}
