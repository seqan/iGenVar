#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "cli_test.hpp"

TEST_F(find_deletions, no_options)
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

TEST_F(find_deletions, fail_no_argument)
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

TEST_F(find_deletions, with_one_argument)
{
    cli_test_result result = execute_app("find_deletions",
                                         "-i", data("detect_breakends_shorted.vcf"));
    std::string expected
    {
        "##fileformat=VCFv4.3\n"
        "##source=iGenVarCaller\n"
        "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
        "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
        "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
        "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr21\t9435236\t.\tN\t<DEL>\t60\tPASS\tEND=9435286;SVLEN=-50;SVTYPE=DEL\n"
        "chr21\t11104574\t.\tN\t<DEL>\t60\tPASS\tEND=11104822;SVLEN=-248;SVTYPE=DEL\n"
    };

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(find_deletions, with_arguments)
{
    cli_test_result result = execute_app("find_deletions",
                                         "-i", data("detect_breakends_shorted.vcf"),
                                         "-o", "find_deletions_file_out_test.vcf");

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, std::string{});
}
