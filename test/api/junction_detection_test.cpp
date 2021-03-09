#include <gtest/gtest.h>

#include "variant_detection/variant_detection.hpp"

std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
std::filesystem::path empty_path{};
const std::vector<detection_methods> default_methods{cigar_string, split_read, read_pairs, read_depth};
const uint64_t sv_default_length = 30;
const uint64_t sample_size = 10;

// Explanation for the strings:
// Reference\tchr21\t41972616\tForward\tRead\t0\t2294\tForward\t1
// INS from Primary Read - Sequence Type: Reference; Chromosome: chr21; Position: 41972616; Orientation: Forward
//                         Sequence Type: Read; Chromosome: 0; Position: 2294; Orientation: Forward
//                         Supporting Reads: 1
// Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\t1
// BND from SA Tag - Sequence Type: Reference; Chromosome: chr22; Position: 17458417; Orientation: Forward
//                   Sequence Type: Reference; Chromosome: chr21; Position: 41972615; Orientation: Forward
//                   Sequence Name: m41327/11677/CCS
// std::string expected_res
// {
//     "Reference\tchr21\t41972616\tForward\tRead\t0\t2294\tForward\t1\n"
//     "Reference\tchr21\t41972616\tReverse\tRead\t0\t3975\tReverse\t1\n"
//     "Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\t1\n"
//     "Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\t2\n"
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

void check_output_and_cleanup(std::string expected_res)
{
    std::string std_cout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(expected_res, std_cout);

    // cleanup
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove output

    return;
}

TEST(junction_detection, fasta_out_not_empty)
{
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       default_methods,
                                       simple_clustering,
                                       no_refinement,
                                       sv_default_length,
                                       empty_path,
                                       sample_size);

    // check_output_and_cleanup(expected_res);
    check_output_and_cleanup(empty_res);
}

/* -------- clustering methods tests -------- */

TEST(junction_detection, clustering_method_hierarchical)
{
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       default_methods,
                                       hierarchical_clustering,
                                       no_refinement,
                                       sv_default_length,
                                       empty_path,
                                       sample_size);

    // check_output_and_cleanup("");
    check_output_and_cleanup(empty_res);
}

TEST(junction_detection, clustering_method_self_balancing_binary_tree)
{
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       default_methods,
                                       self_balancing_binary_tree,
                                       no_refinement,
                                       sv_default_length,
                                       empty_path,
                                       sample_size);

    // check_output_and_cleanup("");
    check_output_and_cleanup(empty_res);
}

TEST(junction_detection, clustering_method_candidate_selection_based_on_voting)
{
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       default_methods,
                                       candidate_selection_based_on_voting,
                                       no_refinement,
                                       sv_default_length,
                                       empty_path,
                                       sample_size);

    // check_output_and_cleanup("");
    check_output_and_cleanup(empty_res);
}

/* -------- refinement methods tests -------- */

TEST(junction_detection, refinement_method_sViper)
{
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       default_methods,
                                       simple_clustering,
                                       sViper_refinement_method,
                                       sv_default_length,
                                       empty_path,
                                       sample_size);

    // check_output_and_cleanup(expected_res);
    check_output_and_cleanup(empty_res);
}

TEST(junction_detection, refinement_method_sVirl)
{
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       default_methods,
                                       simple_clustering,
                                       sVirl_refinement_method,
                                       sv_default_length,
                                       empty_path,
                                       sample_size);

    // check_output_and_cleanup(expected_res);
    check_output_and_cleanup(empty_res);
}
