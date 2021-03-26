#include <fstream>
#include <gtest/gtest.h>

#include <seqan3/io/exception.hpp>

#include "variant_detection/variant_detection.hpp"

const std::string empty_alignment_short_reads_file_path = "";
const std::string default_alignment_long_reads_file_path = DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam";
const std::filesystem::path empty_output_path{};
const std::vector<detection_methods> default_methods{cigar_string, split_read, read_pairs, read_depth};
const uint64_t sv_default_length = 30;

// Explanation for the strings:
// chr21\t41972615\tForward\tchr21\t41972616\tForward\t1\t1681
// INS - Chromosome: chr21; Position: 41972615; Orientation: Forward
//       Chromosome: chr21; Position: 41972616; Orientation: Forward
//       Supporting Reads: 1
//       Average insertion size between two positions: 1681bp
// chr22\t17458417\tForward\tchr21\t41972615\tForward\t1\t2
// BND - Chromosome: chr22; Position: 17458417; Orientation: Forward
//       Chromosome: chr21; Position: 41972615; Orientation: Forward
//       Supporting Reads: 1
//       Average insertion size between two positions: 2bp
// std::string expected_res
// {
//     "chr21\t41972615\tForward\tchr21\t41972616\tForward\t1\t1681\n"
//     "chr22\t17458417\tForward\tchr21\t41972615\tForward\t1\t2\n"
//     "chr22\t17458418\tForward\tchr21\t41972616\tForward\t2\t0\n"
// };

std::string empty_res
{
    "##fileformat=VCFv4.3\n"
    "##source=iGenVarCaller\n"
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
};


void check_output(std::string expected_res)
{
    std::string std_cout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(expected_res, std_cout);

    return;
}

TEST(junction_detection, fasta_out_not_empty)
{
    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                      default_alignment_long_reads_file_path,
                                      default_methods,
                                      simple_clustering,
                                      no_refinement,
                                      sv_default_length,
                                      empty_output_path);

    // check_output_and_cleanup(expected_res);
    check_output(empty_res);
}

TEST(junction_detection, sam_file_unsorted)
{
    // Create a blank SAM file without a sorting indicator.
    const std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::path unsorted_sam_path{tmp_dir/"unsorted.sam"};
    std::ofstream unsorted_sam{unsorted_sam_path.c_str()};
    unsorted_sam << "@HD\tVN:1.6\n" <<
                    "@SQ\tSN:testchr\tLN:1000\n" <<
                    "test1\t16\ttestchr\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n";
    unsorted_sam.close();
    EXPECT_THROW(detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                                   unsorted_sam_path,
                                                   default_methods,
                                                   simple_clustering,
                                                   no_refinement,
                                                   sv_default_length,
                                                   empty_output_path), seqan3::format_error);

    std::filesystem::remove(unsorted_sam_path);
}
/* -------- clustering methods tests -------- */

TEST(junction_detection, clustering_method_hierarchical)
{
    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                      default_alignment_long_reads_file_path,
                                      default_methods,
                                      hierarchical_clustering,
                                      no_refinement,
                                      sv_default_length,
                                      empty_output_path);

    // check_output_and_cleanup("");
    check_output(empty_res);
}

TEST(junction_detection, clustering_method_self_balancing_binary_tree)
{
    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                      default_alignment_long_reads_file_path,
                                      default_methods,
                                      self_balancing_binary_tree,
                                      no_refinement,
                                      sv_default_length,
                                      empty_output_path);

    // check_output_and_cleanup("");
    check_output(empty_res);
}

TEST(junction_detection, clustering_method_candidate_selection_based_on_voting)
{
    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                      default_alignment_long_reads_file_path,
                                      default_methods,
                                      candidate_selection_based_on_voting,
                                      no_refinement,
                                      sv_default_length,
                                      empty_output_path);

    // check_output_and_cleanup("");
    check_output(empty_res);
}

/* -------- refinement methods tests -------- */

TEST(junction_detection, refinement_method_sViper)
{
    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                      default_alignment_long_reads_file_path,
                                      default_methods,
                                      simple_clustering,
                                      sViper_refinement_method,
                                      sv_default_length,
                                      empty_output_path);

    // check_output_and_cleanup(expected_res);
    check_output(empty_res);
}

TEST(junction_detection, refinement_method_sVirl)
{
    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                      default_alignment_long_reads_file_path,
                                      default_methods,
                                      simple_clustering,
                                      sVirl_refinement_method,
                                      sv_default_length,
                                      empty_output_path);

    // check_output_and_cleanup(expected_res);
    check_output(empty_res);
}
