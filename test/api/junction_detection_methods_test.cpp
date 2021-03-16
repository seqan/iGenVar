#include <gtest/gtest.h>

#include <seqan3/core/debug_stream.hpp>

#include "variant_detection/variant_detection.hpp"

std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
std::filesystem::path empty_path{};
const std::vector<detection_methods> all_methods{cigar_string, split_read, read_pairs, read_depth};
const uint64_t sv_default_length = 30;
const uint64_t sample_size = 10;

// Explanation for the strings:
// Reference\tm2257/8161/CCS\t41972616\tForward\tRead\t0\t2294\tForward\tchr21
// INS from Primary Read - Sequence Type: Reference; Sequence Name: m2257/8161/CCS; Position: 41972616; Orientation: Reverse
//                         Sequence Type: Read; Sequence Name: 0; Position: 3975; Orientation: Reverse
//                         Chromosome: chr21
// std::string expected_res_cigar
// {
//     "Reference\tchr21\t41972616\tForward\tRead\t0\t2294\tForward\t1\n"
//     "Reference\tchr21\t41972616\tReverse\tRead\t0\t3975\tReverse\t1\n"
// }

// Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\tm41327/11677/CCS
// BND from SA Tag - Sequence Type: Reference; Chromosome: chr22; Position: 17458417; Orientation: Forward
//                   Sequence Type: Reference; Chromosome: chr21; Position: 41972615; Orientation: Forward
//                   Sequence Name: m41327/11677/CCS
// std::string expected_res_split
// {
//     "Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\t1\n"
//     "Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\t2\n"
// };

// std::string expected_res_pair = "";
// std::string expected_res_depth = "";

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

/* -------- detection methods tests -------- */
/* -------- single methods -------- */

TEST(junction_detection, single_method_only)
{
    for (detection_methods method : all_methods)
    {
        std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");    // remove old output if existent

        testing::internal::CaptureStdout();
        seqan3::debug_stream << "-----------------------------------------------------------------------\n"
                             << "Test Method: " << method << '\n'
                             << "-----------------------------------------------------------------------\n";
        detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                           tmp_dir/"detect_breakends_out_short.fasta",
                                           {method},
                                           simple_clustering,
                                           no_refinement,
                                           sv_default_length,
                                           empty_path,
                                           sample_size);

        //TODO (eldariont): Currently, this (CLI-like) test compares the stdout of the method with an expected string in VCF format.
        //We need to replace this with real API tests that check inidividual parts of the pipeline.
        if (method == 0)
        {
            // check_output_and_cleanup(expected_res_cigar);
            check_output_and_cleanup(empty_res);
        }
        else if (method == 1)
        {
            // check_output_and_cleanup(expected_res_split);
            check_output_and_cleanup(empty_res);
        }
        else if (method == 2)
        {
            // check_output_and_cleanup(expected_res_pair);
            check_output_and_cleanup(empty_res);
        }
        else // (method == 3)
        {
            // check_output_and_cleanup(expected_res_depth);
            check_output_and_cleanup(empty_res);
        }
    }
}

/* -------- method pairs -------- */

TEST(junction_detection, method_pairs)
{
    for (detection_methods method_i : all_methods)
    {
        for (detection_methods method_j : all_methods)
        {
            if (method_i >= method_j) continue; // only check in one order to safe time
            std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");    // remove old output if existent

            testing::internal::CaptureStdout();
            seqan3::debug_stream << "-----------------------------------------------------------------------\n"
                                    << "Test Methods: " << method_i << ", " << method_j << '\n'
                                    << "-----------------------------------------------------------------------\n";
            detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                               tmp_dir/"detect_breakends_out_short.fasta",
                                               {method_i, method_j},
                                               simple_clustering,
                                               no_refinement,
                                               sv_default_length,
                                               empty_path,
                                               sample_size);

            //TODO (eldariont): Currently, this (CLI-like) test compares the stdout of the method with an expected string in VCF format.
            //We need to replace this with real API tests that check inidividual parts of the pipeline.
            if (method_i == 0 && method_j == 1)
            {
                // check_output_and_cleanup(expected_res_cigar + expected_res_split);
                check_output_and_cleanup(empty_res);
            }
            else if (method_i == 0 && method_j == 2)
            {
                // check_output_and_cleanup(expected_res_cigar + expected_res_pair);
                check_output_and_cleanup(empty_res);
            }
            else if (method_i == 0 && method_j == 3)
            {
                // check_output_and_cleanup(expected_res_cigar + expected_res_depth);
                check_output_and_cleanup(empty_res);
            }
            else if (method_i == 1 && method_j == 2)
            {
                // check_output_and_cleanup(expected_res_split + expected_res_pair);
                check_output_and_cleanup(empty_res);
            }
            else if (method_i == 1 && method_j == 3)
            {
                // check_output_and_cleanup(expected_res_split + expected_res_depth);
                check_output_and_cleanup(empty_res);
            }
            else // (method_i == 2 && method_j == 3)
            {
                // check_output_and_cleanup(expected_res_pair + expected_res_depth);
                check_output_and_cleanup(empty_res);
            }
        }
    }
}

/* -------- method triples -------- */

TEST(junction_detection, method_triples)
{
    for (detection_methods method_i : all_methods)
    {
        for (detection_methods method_j : all_methods)
        {
            if (method_i >= method_j) continue; // only check in one order to safe time
            for (detection_methods method_k : all_methods)
            {
                if (method_j >= method_k) continue; // only check in one order to safe time
                std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta"); // remove old output if existent

                testing::internal::CaptureStdout();
                seqan3::debug_stream << "-----------------------------------------------------------------------\n"
                                    << "Test Methods: " << method_i << ", " << method_j << ", " << method_k
                                    << '\n'
                                    << "-----------------------------------------------------------------------\n";
                detect_variants_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                                   tmp_dir/"detect_breakends_out_short.fasta",
                                                   {method_i, method_j, method_k},
                                                   simple_clustering,
                                                   no_refinement,
                                                   sv_default_length,
                                                   empty_path,
                                                   sample_size);

                //TODO (eldariont): Currently, this (CLI-like) test compares the stdout of the method with an expected string in VCF format.
                //We need to replace this with real API tests that check inidividual parts of the pipeline.
                if (method_i == 0 && method_j == 1 && method_k == 2)
                {
                    // check_output_and_cleanup(expected_res_cigar + expected_res_split + expected_res_pair);
                    check_output_and_cleanup(empty_res);
                }
                else if (method_i == 0 && method_j == 1 && method_k == 3)
                {
                    // check_output_and_cleanup(expected_res_cigar + expected_res_split + expected_res_depth);
                    check_output_and_cleanup(empty_res);
                }
                else if (method_i == 0 && method_j == 2 && method_k == 3)
                {
                    // check_output_and_cleanup(expected_res_cigar + expected_res_pair + expected_res_depth);
                    check_output_and_cleanup(empty_res);
                }
                else if (method_i == 1 && method_j == 2 && method_k == 3)
                {
                    // check_output_and_cleanup(expected_res_split + expected_res_pair + expected_res_depth);
                    check_output_and_cleanup(empty_res);
                }
            }
        }
    }
}
