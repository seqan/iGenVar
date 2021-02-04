#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan3/core/debug_stream.hpp>

#include "detect_breakends/junction_detection.hpp"

// TEST(group1, fasta_out_empty)
// {
//     std::string expected{"Reference\tchr9\t70103073\tForward\tReference\tchr9\t70103147\tForward\tm13802/6999/CCS\t1\n"};
//     testing::internal::CaptureStdout();
//     detect_junctions_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam", "");
//     std::string std_cout = testing::internal::GetCapturedStdout();
//     EXPECT_EQ(expected, std_cout);
// }

// Explanation for the stings:
// Reference\tm2257/8161/CCS\t41972616\tForward\tRead\t0\t2294\tForward\tchr21
// INS from Primary Read - Sequence Type: Reference; Sequence Name: m2257/8161/CCS; Position: 41972616; Orientation: Reverse
//                         Sequence Type: Read; Sequence Name: 0; Position: 3975; Orientation: Reverse
//                         Chromosome: chr21
// Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\tm41327/11677/CCS
// BND from SA Tag - Sequence Type: Reference; Chromosome: chr22; Position: 17458417; Orientation: Forward
//                   Sequence Type: Reference; Chromosome: chr21; Position: 41972615; Orientation: Forward
//                   Sequence Name: m41327/11677/CCS

uint64_t sv_default_length = 30;

TEST(junction_detection, fasta_out_not_empty)
{
    std::string expected{
        "Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\t1\n"
        "Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\t2\n"
        "Reference\tm2257/8161/CCS\t41972616\tForward\tRead\t0\t2294\tForward\t1\n"
        "Reference\tm2257/8161/CCS\t41972616\tReverse\tRead\t0\t3975\tReverse\t1\n"
    };

    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_junctions_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       {1, 2, 3, 4},
                                       simple_clustering,
                                       no_refinement,
                                       sv_default_length);

    std::string std_cout = testing::internal::GetCapturedStdout();
    seqan3::debug_stream << "std_out:\n" << std_cout << '\n';
    EXPECT_EQ(expected, std_cout);

    // cleanup
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove output
}

TEST(junction_detection, method_1_only)
{
    std::string expected{
        "Reference\tm2257/8161/CCS\t41972616\tForward\tRead\t0\t2294\tForward\t1\n"
        "Reference\tm2257/8161/CCS\t41972616\tReverse\tRead\t0\t3975\tReverse\t1\n"
    };

    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_junctions_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       {1},
                                       simple_clustering,
                                       no_refinement,
                                       sv_default_length);

    std::string std_cout = testing::internal::GetCapturedStdout();
    seqan3::debug_stream << "std_out:\n" << std_cout << '\n';
    EXPECT_EQ(expected, std_cout);

    // cleanup
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove output
}

TEST(junction_detection, method_2_only)
{
    std::string expected{
        "Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\t1\n"
        "Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\t2\n"
    };

    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_junctions_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       {2},
                                       simple_clustering,
                                       no_refinement,
                                       sv_default_length);

    std::string std_cout = testing::internal::GetCapturedStdout();
    seqan3::debug_stream << "std_out:\n" << std_cout << '\n';
    EXPECT_EQ(expected, std_cout);

    // cleanup
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove output
}

TEST(junction_detection, method_1_and_2)
{
    std::string expected{
        "Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\t1\n"
        "Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\t2\n"
        "Reference\tm2257/8161/CCS\t41972616\tForward\tRead\t0\t2294\tForward\t1\n"
        "Reference\tm2257/8161/CCS\t41972616\tReverse\tRead\t0\t3975\tReverse\t1\n"
    };

    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_junctions_in_alignment_file(DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam",
                                       tmp_dir/"detect_breakends_out_short.fasta",
                                       {1, 2},
                                       simple_clustering,
                                       no_refinement,
                                       sv_default_length);

    std::string std_cout = testing::internal::GetCapturedStdout();
    seqan3::debug_stream << "std_out:\n" << std_cout << '\n';
    EXPECT_EQ(expected, std_cout);

    // cleanup
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove output
}
