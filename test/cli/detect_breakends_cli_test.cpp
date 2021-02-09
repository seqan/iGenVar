#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "cli_test.hpp"
#include <fstream>
#include <sstream>
TEST_F(detect_breakends, no_options)
{
    cli_test_result result = execute_app("detect_breakends");
    std::string expected
    {
            "iGenVar - Detect junctions in a read alignment file\n"
            "===================================================\n"
            "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(detect_breakends, fail_no_argument)
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


TEST_F(detect_breakends, with_arguments)
{
    cli_test_result result = execute_app("detect_breakends",
                                         data("simulated.minimap2.hg19.coordsorted_cutoff.sam"),
                                         "detect_breakends_insertion_file_out.fasta");
    std::string expected
    {
        "Reference\tchr21\t41972616\tForward\tRead\t0\t2294\tForward\t1\n"
        "Reference\tchr21\t41972616\tReverse\tRead\t0\t3975\tReverse\t1\n"
        "Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\t1\n"
        "Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\t2\n"
    };
    std::string expected_err
    {
        "INS1: Reference\tchr21\t41972616\tForward\tRead\t0\t2294\tForward\tm2257/8161/CCS\n"
        "INS2: Reference\tchr21\t41972616\tReverse\tRead\t0\t3975\tReverse\tm2257/8161/CCS\n"
        "The read pair method is not yet implemented.\n"
        "The read depth method is not yet implemented.\n"
        "BND: Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\tm41327/11677/CCS\n"
        "The read pair method is not yet implemented.\n"
        "The read depth method is not yet implemented.\n"
        "BND: Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\tm21263/13017/CCS\n"
        "The read pair method is not yet implemented.\n"
        "The read depth method is not yet implemented.\n"
        "BND: Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\tm38637/7161/CCS\n"
        "The read pair method is not yet implemented.\n"
        "The read depth method is not yet implemented.\n"
        "Start clustering...\n"
        "Done with clustering. Found 4 junction clusters.\n"
        "No refinement was selected.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(detect_breakends, test_outfile)
{
    cli_test_result result = execute_app("detect_breakends",
                                         data("simulated.minimap2.hg19.coordsorted_cutoff.sam"),
                                         "detect_breakends_insertion_file_out.fasta");

    std::cout << "Current path is " << std::filesystem::current_path() << '\n';

    std::filesystem::path out_file_path = "detect_breakends_insertion_file_out.fasta";
    std::filesystem::path test_file_path = "../../data/detect_breakends_insertion_file_test.fasta";
    seqan3::sequence_file_input out_file{out_file_path};
    seqan3::sequence_file_input test_file{test_file_path};

    EXPECT_RANGE_EQ(out_file, test_file);
}
