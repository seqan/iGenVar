#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "cli_test.hpp"
#include <fstream>
#include <sstream>
TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("detect_breakends");
    std::string expected
    {
            "detectJunctions - Detect junctions in a read alignment file\n"
            "===========================================================\n"
            "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("detect_breakends", "-v");
    std::string expected
    {
        "[Error] Unknown option -v. In case this is meant to be a non-option/argument/parameter, please specify "
        "the start of non-options with '--'. See -h/--help for program information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}


TEST_F(cli_test, with_arguments)
{
    cli_test_result result = execute_app("detect_breakends",
                                         data("converted_bam_shorted.sam"),
                                         "detect_breakends_insertion_file_out.fasta");
    std::string expected
    {
        "Reference\tchr9\t70103073\tForward\tReference\tchr9\t70103147\tForward\tm13802/6999/CCS\n"
    };
    std::string expected_err
    {
        "DEL: Reference\tchr9\t70103073\tForward\tReference\tchr9\t70103147\tForward\tm13802/6999/CCS\nDone. Found 1 junctions.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(cli_test, test_outfile)
{
    cli_test_result result = execute_app("detect_breakends",
                                         data("converted_bam_shorted.sam"),
                                         "detect_breakends_insertion_file_out.fasta");
    std::ifstream f;
    f.open("detect_breakends_insertion_file_out.fasta");
    std::stringstream buffer;
    buffer << f.rdbuf();
    //expected string currently empty:
    std::string expected
    {
	    ""
    };
    //this does not specifically check if file exists, rather if its readable.
    EXPECT_TRUE(f.is_open());
    EXPECT_EQ(buffer.str(), expected);
}
