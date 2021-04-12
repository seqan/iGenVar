#pragma once

#include <seqan3/argument_parser/argument_parser.hpp>   // for seqan3::argument_parser

#include "variant_detection/variant_detection.hpp"

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args);

/*! \brief Detects genomic variants by analyzing an alignment file (sam/bam). The detected
 *         variants are printed to a given file or stdout and insertion alleles are stored in a fasta file.
 *
 * \param[in] args - command line arguments:\n
 *                   **args.alignment_short_reads_file_path** - short reads input file, path to the sam/bam file\n
 *                   **args.alignment_long_reads_file_path** - long reads input file, path to the sam/bam file\n
 *                   **args.output_file_path** output file - path for the VCF file - *default: standard output*\n
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
 *                   **args.min_var_length** - minimum length of variants to detect - *default: 30 bp*
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
