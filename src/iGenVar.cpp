#include "iGenVar.hpp"

#include <map>

#include <seqan3/contrib/stream/bgzf_stream_util.hpp>       // for bgzf_thread_count
#include <seqan3/core/debug_stream.hpp>                     // for seqan3::debug_stream
#include <seqan3/io/sequence_file/format_embl.hpp>          // for embl sequence file extensions
#include <seqan3/io/sequence_file/format_fasta.hpp>         // for fasta sequence file extensions
#include <seqan3/io/sequence_file/format_fastq.hpp>         // for fastq sequence file extensions
#include <seqan3/io/sequence_file/format_genbank.hpp>       // for genbank sequence file extensions

#include "modules/clustering/hierarchical_clustering_method.hpp"    // for the hierarchical clustering method
#include "modules/clustering/simple_clustering_method.hpp"          // for the simple clustering method
#include "structures/cluster.hpp"                                   // for class Cluster
#include "variant_detection/snp_indel_detection.hpp"                // for detect_snp_and_indel
#include "variant_detection/variant_detection.hpp"                  // for detect_junctions_in_long_reads_sam_file()
#include "variant_detection/variant_output.hpp"                     // for find_and_output_variants()

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "Lydia Buntrock, David Heller, Joshua Kim";
    parser.info.app_name = "iGenVar";
    parser.info.man_page_title = "Short and Long Read SV Caller";
    parser.info.short_description = "Detect genomic variants in a read alignment file";
    parser.info.version = "0.0.3";
    parser.info.date = "30-03-2021";    // last update
    parser.info.email = "lydia.buntrock@fu-berlin.de";
    parser.info.long_copyright = "long_copyright";
    parser.info.short_copyright = "short_copyright";
    parser.info.url = "https://github.com/seqan/iGenVar/";

    // merge the vectors of all sequence file extensions
    std::vector<std::string> seq_file_extensions = seqan3::format_fasta::file_extensions;
    seq_file_extensions.insert(seq_file_extensions.end(),
                               seqan3::format_fastq::file_extensions.begin(),
                               seqan3::format_fastq::file_extensions.end());
    seq_file_extensions.insert(seq_file_extensions.end(),
                               seqan3::format_embl::file_extensions.begin(),
                               seqan3::format_embl::file_extensions.end());
    seq_file_extensions.insert(seq_file_extensions.end(),
                               seqan3::format_genbank::file_extensions.begin(),
                               seqan3::format_genbank::file_extensions.end());

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
    parser.add_option(args.genome_file_path,
                      'g', "input_genome",
                      "Input the reference genome in FASTA or FASTQ format.",
                      seqan3::option_spec::standard,
                      seqan3::input_file_validator{seq_file_extensions} );
    parser.add_option(args.output_file_path, 'o', "output",
                      "The path of the vcf output file. If no path is given, will output to standard output.",
                      seqan3::option_spec::standard,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create, {"vcf"}});

    // Options - VCF:
    parser.add_option(args.vcf_sample_name, 's', "vcf_sample_name",
                      "Specify your sample name for the vcf header line.",
                      seqan3::option_spec::standard);

    // Options - Other parameters:
    parser.add_option(args.threads, 't', "threads",
                      "Specify the number of decompression threads used for reading BAM files.",
                      seqan3::option_spec::standard);
    parser.add_flag(gVerbose, 'v', "verbose",
                    "If you set this flag, we provide additional details about what iGenVar does. The detailed output "
                    "is printed in the standard error.");

    // Options - Optional output:
    parser.add_option(args.junctions_file_path, 'a', "junctions",
                      "The path of the optional junction output file. If no path is given, junctions will not be output.",
                      seqan3::option_spec::advanced,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create});
    parser.add_option(args.clusters_file_path, 'b', "clusters",
                      "The path of the optional cluster output file. If no path is given, clusters will not be output.",
                      seqan3::option_spec::advanced,
                      seqan3::output_file_validator{seqan3::output_file_open_options::open_or_create});

    // Options - Methods:
    parser.add_option(args.methods, 'd', "method",
                      "Choose the detection method(s) to be used. "
                      "Value must be one of (method name or number) "
                      "[0,cigar_string,1,split_read,2,read_pairs,3,read_depth].",
                      seqan3::option_spec::advanced);
    parser.add_option(args.clustering_method, 'c', "clustering_method",
                      "Choose the clustering method to be used. "
                      "Value must be one of (method name or number) [0,simple_clustering,"
                      "1,hierarchical_clustering,2,self_balancing_binary_tree,3,candidate_selection_based_on_voting].",
                      seqan3::option_spec::advanced);
    parser.add_option(args.refinement_method, 'r', "refinement_method",
                      "Choose the refinement method to be used. "
                      "Value must be one of (method name or number) "
                      "[0,no_refinement,1,sViper_refinement_method,2,sVirl_refinement_method].",
                      seqan3::option_spec::advanced);

    // Options - SV specifications:
    parser.add_option(args.min_var_length, 'k', "min_var_length",
                      "Specify what should be the minimum length of your SVs to be detected. "
                      "This value needs to be non-negative.",
                      seqan3::option_spec::advanced);
    parser.add_option(args.max_var_length, 'l', "max_var_length",
                      "Specify what should be the maximum length of your SVs to be detected. "
                      "SVs larger than this threshold can still be output as translocations. "
                      "This value needs to be non-negative.",
                      seqan3::option_spec::advanced);
    parser.add_option(args.max_tol_inserted_length, 'm', "max_tol_inserted_length",
                      "Specify what should be the longest tolerated inserted sequence at sites of non-INS SVs. "
                      "This value needs to be non-negative.",
                      seqan3::option_spec::advanced);
    parser.add_option(args.max_tol_deleted_length, 'e', "max_tol_deleted_length",
                      "Specify what should be the longest tolerated deleted sequence at sites of non-DEL SVs. "
                      "This value needs to be non-negative.",
                      seqan3::option_spec::advanced);
    parser.add_option(args.max_overlap, 'n', "max_overlap",
                      "Specify the maximum allowed overlap between two alignment segments. "
                      "This value needs to be non-negative.",
                      seqan3::option_spec::advanced);
    parser.add_option(args.min_qual, 'q', "min_qual",
                      "Specify the minimum quality (amount of supporting reads) of a structural variant to be reported "
                      "in the vcf output file. This value needs to be non-negative.",
                      seqan3::option_spec::advanced);

    // Options - Clustering specifications:
    parser.add_option(args.partition_max_distance, 'p', "partition_max_distance",
                      "Specify the maximum distance in bp between members of the same partition."
                      "This value needs to be non-negative.",
                      seqan3::option_spec::advanced);
    parser.add_option(args.hierarchical_clustering_cutoff, 'w', "hierarchical_clustering_cutoff",
                      "Specify the distance cutoff for the hierarchical clustering. "
                      "This value needs to be non-negative.",
                      seqan3::option_spec::advanced);
}

void detect_variants_in_alignment_file(cmd_arguments const & args)
{
    // Store junctions
    std::vector<Junction> junctions{};
    std::map<std::string, int32_t> references_lengths{};

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

    if (args.junctions_file_path != "")
    {
        std::ofstream junctions_file{args.junctions_file_path};
        if (!junctions_file.good() || !junctions_file.is_open())
        {
            throw std::runtime_error{"Could not open file '" + args.junctions_file_path.string() + "' for writing."};
        }
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

    if (args.clusters_file_path != "")
    {
        std::ofstream clusters_file{args.clusters_file_path};
        if (!clusters_file.good() || !clusters_file.is_open())
        {
            throw std::runtime_error{"Could not open file '" + args.clusters_file_path.string() + "' for writing."};
        }
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
