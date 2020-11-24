#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "cli_test.hpp"

TEST_F(cli_test, no_options)
{
    cli_test_result result = execute_app("find_deletions");
    std::string expected
    {
        "partitionJunctions - Find deletions\n"
        "===================================\n"
        "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(cli_test, fail_no_argument)
{
    cli_test_result result = execute_app("find_deletions", "-v");
    std::string expected
    {
        "[Error] Option -i/--input is required but not set.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(cli_test, with_one_argument)
{
    cli_test_result result = execute_app("find_deletions",
                                         "-i", data("detect_breakends_shorted.vcf"));
    
    std::string err
    {
        "[Error] Option -o/--output is required but not set.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, err);
}

TEST_F(cli_test, with_arguments)
{
    cli_test_result result = execute_app("find_deletions",
                                         "-i", data("detect_breakends_shorted.vcf"),
                                         "-o", "find_deletions_file_out_test.vcf");
    
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}
