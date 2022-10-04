#include "iGenVar.hpp"

#include <map>

#include <seqan3/contrib/stream/bgzf_stream_util.hpp>       // for bgzf_thread_count
#include <seqan3/core/debug_stream.hpp>                     // for seqan3::debug_stream
#include <seqan3/io/sequence_file/input.hpp>                // for sequence input file extensions

#include "modules/clustering/hierarchical_clustering_method.hpp"    // for the hierarchical clustering method
#include "modules/clustering/simple_clustering_method.hpp"          // for the simple clustering method
#include "structures/cluster.hpp"                                   // for class Cluster
#include "variant_detection/snp_indel_detection.hpp"                // for detect_snp_and_indel
#include "variant_detection/variant_detection.hpp"                  // for detect_junctions_in_long_reads_sam_file()
#include "variant_detection/variant_output.hpp"                     // for find_and_output_variants()

static inline std::vector<std::string> sequence_extensions{
    seqan3::detail::valid_file_extensions<typename seqan3::sequence_file_input<>::valid_formats>()};

static inline std::vector<std::string> compression_extensions{[]()
                                                              {
                                                                  std::vector<std::string> result;
#ifdef SEQAN3_HAS_BZIP2
                                                                  result.push_back("bz2");
#endif
#ifdef SEQAN3_HAS_ZLIB
                                                                  result.push_back("gz");
                                                                  result.push_back("bgzf");
#endif
                                                                  return result;
                                                              }()}; // GCOVR_EXCL_LINE

static inline std::vector<std::string> combined_extensions{
    []()
    {
        if (compression_extensions.empty())
            return sequence_extensions; // GCOVR_EXCL_LINE
        std::vector<std::string> result;
        for (auto && sequence_extension : sequence_extensions)
        {
            result.push_back(sequence_extension);
            for (auto && compression_extension : compression_extensions)
                result.push_back(sequence_extension + std::string{'.'} + compression_extension);
        }
        return result;
    }()};

void initialize_argument_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Lydia Buntrock, David Heller, Joshua Kim";
    parser.info.app_name = "iGenVar";
    parser.info.man_page_title = "Short and Long Read SV Caller";
    parser.info.short_description = "Detect genomic variants in a read alignment file";
    parser.info.version = "0.0.3";
    parser.info.date = "11-07-2022";    // last update
    parser.info.email = "lydia.buntrock@fu-berlin.de";
    parser.info.long_copyright = "long_copyright";
    parser.info.short_copyright = "short_copyright";
    parser.info.url = "https://github.com/seqan/iGenVar/";

    // Options - Input / Output:
    parser.add_subsection("Input / Output");
    parser.add_option(args.alignment_short_reads_file_path, sharg::config{
                        .short_id = 'i', .long_id = "input_short_reads",
                        .description = "Input short read alignments in SAM or BAM format (Illumina).",
                        .validator = sharg::input_file_validator{{"sam", "bam"}} });
    parser.add_option(args.alignment_long_reads_file_path, sharg::config{
                        .short_id = 'j', .long_id = "input_long_reads",
                        .description = "Input long read alignments in SAM or BAM format (PacBio, Oxford Nanopore, ...).",
                        .validator = sharg::input_file_validator{{"sam", "bam"}} });
    parser.add_option(args.genome_file_path, sharg::config{
                        .short_id = 'g', .long_id = "input_genome",
                        .description = "Input the sequence of the reference genome.",
                        .validator = sharg::input_file_validator{combined_extensions}});
    parser.add_option(args.output_file_path, sharg::config{
                        .short_id = 'o', .long_id = "output",
                        .description = "The path of the vcf output file. If no path is given, will output to standard output.",
                        .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"vcf"}}});

    // Options - VCF:
    parser.add_subsection("VCF");
    parser.add_option(args.vcf_sample_name, sharg::config{
                        .short_id = 's', .long_id = "vcf_sample_name",
                        .description = "Specify your sample name for the vcf header line."});

    // Options - Other parameters:
    parser.add_subsection("Other parameters");
    parser.add_option(args.threads, sharg::config{
                        .short_id = 't', .long_id = "threads",
                        .description = "Specify the number of decompression threads used for reading BAM files."});
    parser.add_flag(gVerbose, sharg::config{
                        .short_id = 'v', .long_id = "verbose",
                        .description = "If you set this flag, we provide additional details about what iGenVar does. "
                                       "The detailed output is printed in the standard error."});

    // Options - Optional output:
    parser.add_subsection("Optional output");
    parser.add_option(args.junctions_file_path, sharg::config{
                        .short_id = 'a', .long_id = "junctions",
                        .description = "The path of the optional junction output file. If no path is given, junctions will not be output.",
                        .advanced = true,
                        .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create}});
    parser.add_option(args.clusters_file_path, sharg::config{
                        .short_id = 'b', .long_id = "clusters",
                        .description = "The path of the optional cluster output file. If no path is given, clusters will not be output.",
                        .advanced = true,
                        .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create}});

    // Options - Methods:
    parser.add_subsection("Methods");
    parser.add_option(args.methods, sharg::config{
                        .short_id = 'd', .long_id = "method",
                        .description = "Choose the detection method(s) to be used. "
                                       "Value must be one of (method name or number) "
                                       "[0,cigar_string,1,split_read,2,read_pairs,3,read_depth].",
                        .advanced = true});
    parser.add_option(args.clustering_method, sharg::config{
                        .short_id = 'c', .long_id = "clustering_method",
                        .description = "Choose the clustering method to be used. Value must be one of "
                                       "(method name or number) [0,simple_clustering,1,hierarchical_clustering,"
                                       "2,self_balancing_binary_tree,3,candidate_selection_based_on_voting].",
                        .advanced = true});
    parser.add_option(args.refinement_method, sharg::config{
                        .short_id = 'r', .long_id = "refinement_method",
                        .description = "Choose the refinement method to be used. "
                                       "Value must be one of (method name or number) "
                                       "[0,no_refinement,1,sViper_refinement_method,2,sVirl_refinement_method].",
                        .advanced = true});

    // Options - SV specifications:
    parser.add_subsection("SV specifications");
    parser.add_option(args.min_var_length, sharg::config{
                        .short_id = 'k', .long_id = "min_var_length",
                        .description = "Specify what should be the minimum length of your SVs to be detected. "
                                       "This value needs to be non-negative.",
                        .advanced = true});
    parser.add_option(args.max_var_length, sharg::config{
                        .short_id = 'l', .long_id = "max_var_length",
                        .description = "Specify what should be the maximum length of your SVs to be detected. "
                                       "SVs larger than this threshold can still be output as translocations. "
                                       "This value needs to be non-negative.",
                        .advanced = true});
    parser.add_option(args.max_tol_inserted_length, sharg::config{
                        .short_id = 'm', .long_id = "max_tol_inserted_length",
                        .description = "Specify what should be the longest tolerated inserted sequence at sites of DEL "
                                       "SVs. This value needs to be non-negative.",
                        .advanced = true});
    parser.add_option(args.max_tol_deleted_length, sharg::config{
                        .short_id = 'e', .long_id = "max_tol_deleted_length",
                        .description = "Specify what should be the longest tolerated deleted sequence at sites of "
                                       "non-DEL/INV SVs. This value needs to be non-negative.",
                        .advanced = true});
    parser.add_option(args.max_overlap, sharg::config{
                        .short_id = 'n', .long_id = "max_overlap",
                        .description = "Specify the maximum allowed overlap between two alignment segments. "
                                       "This value needs to be non-negative.",
                        .advanced = true});
    parser.add_option(args.min_qual, sharg::config{
                        .short_id = 'q', .long_id = "min_qual",
                        .description = "Specify the minimum quality (amount of supporting reads) of a structural variant "
                                       "to be reported in the vcf output file. This value needs to be non-negative.",
                        .advanced = true});

    // Options - Clustering specifications:
    parser.add_subsection("Clustering specifications");
    parser.add_option(args.partition_max_distance, sharg::config{
                        .short_id = 'p', .long_id = "partition_max_distance",
                        .description = "Specify the maximum distance in bp between members of the same partition."
                                       "This value needs to be non-negative.",
                        .advanced = true});
    parser.add_option(args.hierarchical_clustering_cutoff, sharg::config{
                        .short_id = 'w', .long_id = "hierarchical_clustering_cutoff",
                        .description = "Specify the distance cutoff for the hierarchical clustering. "
                                       "This value needs to be non-negative.",
                        .advanced = true});
}

void detect_variants_in_alignment_file(cmd_arguments const & args)
{
    // Store junctions
    std::vector<Junction> junctions{};
    // Map of contig names and their length (SN and LN tag of @SQ)
    std::map<std::string, int32_t> references_lengths{};

    // if short and long reads are given warn the user about the necessary equality of references
    if (!args.alignment_short_reads_file_path.empty() && !args.alignment_long_reads_file_path.empty())
    {
        seqan3::debug_stream << "You have specified two input files for short and long read data. Note that they should"
                                " be mapped to the same reference, e.g. GRCh37 (hg19) or GRCh38 (hg38). If they come"
                                " from different versions, for example the coordinates may not match. In such a case,"
                                " use a coordinate converter beforehand.\n";
    }

    // short reads
    // TODO (joergi-w 30.09.2021) Control the selection with the 'method' parameter, not the availability of a genome.
    if (!args.alignment_short_reads_file_path.empty() && args.genome_file_path.empty())
    {
        seqan3::debug_stream << "Detect junctions in short reads...\n";
        detect_junctions_in_short_reads_sam_file(junctions, references_lengths, args);
    }

    // long reads
    if (!args.alignment_long_reads_file_path.empty())
    {
        seqan3::debug_stream << "Detect junctions in long reads...\n";
        detect_junctions_in_long_reads_sam_file(junctions, references_lengths, args);
    }

    // SNPs and indels for short reads; genome must be given
    if (!args.alignment_short_reads_file_path.empty() && !args.genome_file_path.empty())
    {
        seqan3::debug_stream << "Detect SNPs and indels in short reads...\n";
        detect_snp_and_indel(args.alignment_short_reads_file_path, args.min_var_length);
    }

    std::sort(junctions.begin(), junctions.end());

    if (!args.junctions_file_path.empty())
    {
        std::ofstream junctions_file{args.junctions_file_path};

        // LCOV_EXCL_START
        if (!junctions_file.good() || !junctions_file.is_open())
            throw std::runtime_error{"Could not open file '" + args.junctions_file_path.string() + "' for writing."};
        // LCOV_EXCL_STOP

        for (Junction const & junction : junctions)
        {
            junctions_file << junction << "\n";
        }
        junctions_file.close();
    }

    seqan3::debug_stream << "Start clustering...\n";

    std::vector<Cluster> clusters;
    switch (args.clustering_method)
    {
        case 0: // simple_clustering
            clusters = simple_clustering_method(junctions);
            break;
        case 1: // hierarchical clustering
            clusters = hierarchical_clustering_method(junctions,
                                                      args.partition_max_distance,
                                                      args.hierarchical_clustering_cutoff);
            break;
        case 2: // self-balancing_binary_tree,
            seqan3::debug_stream << "The self-balancing binary tree clustering method is not yet implemented.\n";
            break;
        case 3: // candidate_selection_based_on_voting
            seqan3::debug_stream << "The candidate selection based on voting clustering method is not yet implemented.\n";
            break;
    }

    seqan3::debug_stream << "Done with clustering. Found " << clusters.size() << " junction clusters.\n";

    if (!args.clusters_file_path.empty())
    {
        std::ofstream clusters_file{args.clusters_file_path};

        // LCOV_EXCL_START
        if (!clusters_file.good() || !clusters_file.is_open())
            throw std::runtime_error{"Could not open file '" + args.clusters_file_path.string() + "' for writing."};
        // LCOV_EXCL_STOP

        for (Cluster const & cluster : clusters)
        {
            clusters_file << cluster << "\n";
        }
        clusters_file.close();
    }

    switch (args.refinement_method)
    {
        case 0: // no refinement
            seqan3::debug_stream << "No refinement was selected.\n";
            break;
        case 1: // sViper_refinement_method
            seqan3::debug_stream << "The sViper refinement method is not yet implemented.\n";
            break;
        case 2: // sVirl_refinement_method
            seqan3::debug_stream << "The sVirl refinement method is not yet implemented.\n";
            break;
    }

    find_and_output_variants(references_lengths, clusters, args, args.output_file_path);
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"iGenVar", argc, argv};                    // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(parser, args);

    // Parse the given arguments and catch possible errors.
    try
    {
        parser.parse();                                             // trigger command line parsing
    }
    catch (sharg::parser_error const & ext)                         // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << '\n';   // customise your error message
        return -1;
    }

    // Check if we have at least one input file.
    if (args.alignment_short_reads_file_path.empty() && args.alignment_long_reads_file_path.empty())
    {
        seqan3::debug_stream << "[Error] You need to input at least one sam/bam file.\n"
                             << "Please use -i or -input_short_reads to pass a short read file "
                                "or -j or -input_long_reads for a long read file.\n";
        return -1;
    }

    // Set the number of decompression threads
    seqan3::contrib::bgzf_thread_count = args.threads;

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

    // Check that the given parameters are non-negative.
    if (args.hierarchical_clustering_cutoff < 0)
    {
        seqan3::debug_stream << "[Error] You gave a negative hierarchical_clustering_cutoff parameter.\n";
        return -1;
    }

    detect_variants_in_alignment_file(args);

    return 0;
}
