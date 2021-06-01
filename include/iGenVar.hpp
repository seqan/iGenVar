#pragma once

#include <seqan3/argument_parser/argument_parser.hpp>   // for seqan3::argument_parser

#include "variant_detection/method_enums.hpp"           // for enum detection_methods, clustering_methods and refinement_methods

struct cmd_arguments
{
    std::filesystem::path alignment_short_reads_file_path{""};
    std::filesystem::path alignment_long_reads_file_path{""};
    std::filesystem::path output_file_path{};
    std::string vcf_sample_name{"GENOTYPE"};
    std::vector<detection_methods> methods{cigar_string, split_read, read_pairs, read_depth};   // default: all methods
    clustering_methods clustering_method{hierarchical_clustering};                              // default: hierarchical clustering method
    refinement_methods refinement_method{no_refinement};                                        // default: no refinement
    int32_t min_var_length = 30;
    int32_t max_var_length = 1000000;
    int32_t max_tol_inserted_length = 5;
    int32_t max_overlap = 10;
    int16_t threads = 1;
    int32_t min_qual = 1;
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args);

/*! \brief Detects genomic variants by analyzing an alignment file (sam/bam). The detected
 *         variants are printed to a given file or stdout and insertion alleles are stored in a fasta file.
 *
 * \param[in] args - command line arguments:\n
 *                   **args.alignment_short_reads_file_path** - short reads input file, path to the sam/bam file\n
 *                   **args.alignment_long_reads_file_path** - long reads input file, path to the sam/bam file\n
 *                   **args.output_file_path** output file - path for the VCF file - *default: standard output*\n
 *                   **args.vcf_sample_name - Name of the sample for the vcf header line*\n
 *                   **args.methods** - list of methods for detecting junctions
 *                      (1: cigar_string, 2: split_read, 3: read_pairs, 4: read_depth) - *default: all methods*\n
 *                   **args.clustering_method** - method for clustering junctions
 *                      (0: simple_clustering,
 *                       1: hierarchical_clustering,
 *                       2: self-balancing_binary_tree,
 *                       3: candidate_selection_based_on_voting) - *default: simple clustering method*\n
 *                   **args.refinement_method** - method for refining junctions
 *                      (0: no_refinement,
 *                       1: sViper_refinement_method,
 *                       2: sVirl_refinement_method) - *default: no refinement*\n
 *                   **args.min_var_length** - minimum length of variants to detect
 *                                             (expected to be non-negative) - *default: 30 bp*\n
 *                   **args.max_var_length** - maximum length of variants to detect
 *                                             (expected to be non-negative) - *default: 1,000,000 bp*\n
 *                   **args.max_tol_inserted_length** - longest tolerated inserted sequence at non-INS SV types
 *                                                      (expected to be non-negative) - *default: 5 bp*\n
 *                   **args.max_overlap** - maximum overlap between alignment segments
 *                                          (expected to be non-negative) - *default: 10 bp*\n
 *                   **args.min_qual** - minimum quality (amount of supporting reads) of a structural variant
 *                                       (expected to be non-negative) - *default: 1 supporting read*\n
 *
 *
 * \details Detects novel junctions from read alignment records using different detection methods.
 *          The junctions are clustered using one of several clustering methods.
 *          Then, the junction clusters are refined using one of several refinement methods.
 *          Finally, the refined junction clusters are categorized into different variant classes
 *          and output in VCF format.
 */
void detect_variants_in_alignment_file(cmd_arguments const & args);

int main(int argc, char ** argv);
