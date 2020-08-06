#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <seqan3/core/debug_stream.hpp>

#include "detect_breakends/junction_detection.hpp"

// TEST(group1, fasta_out_empty)
// {
//     std::string expected{"Reference\tchr9\t70103073\tForward\tReference\tchr9\t70103147\tForward\tm13802/6999/CCS\n"};
//     testing::internal::CaptureStdout();
//     detect_junctions_in_alignment_file(DATADIR"converted_bam_shorted.sam", "");
//     std::string std_cout = testing::internal::GetCapturedStdout();
//     EXPECT_EQ(expected, std_cout);
// }

TEST(junction_detection, fasta_out_not_empty)
{
    std::string expected{"Reference\tchr9\t70103073\tForward\tReference\tchr9\t70103147\tForward\tm13802/6999/CCS\n"};

    std::filesystem::path tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove old output if existent

    testing::internal::CaptureStdout();
    detect_junctions_in_alignment_file(DATADIR"converted_bam_shorted.sam", tmp_dir/"detect_breakends_out_short.fasta");

    std::string std_cout = testing::internal::GetCapturedStdout();
    seqan3::debug_stream << "std_out:\n" << std_cout << '\n';
    EXPECT_EQ(expected, std_cout);

    // cleanup
    std::filesystem::remove(tmp_dir/"detect_breakends_out_short.fasta");        // remove output
}
