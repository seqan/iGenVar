#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>

#include "cli_test.hpp"

std::string help_page
{
    "partitionJunctions - Find deletions\n"
    "===================================\n"
    "\nOPTIONS\n"
    "\n"
    "  Basic options:\n"
    "    -h, --help\n"
    "          Prints the help page.\n"
    "    -hh, --advanced-help\n"
    "          Prints the help page including advanced options.\n"
    "    --version\n"
    "          Prints the version information.\n"
    "    --copyright\n"
    "          Prints the copyright/license information.\n"
    "    --export-help (std::string)\n"
    "          Export the help page information. Value must be one of [html, man].\n"
    "    --version-check (bool)\n"
    "          Whether to check for the newest app version. Default: true.\n"
    "    -i, --input (std::filesystem::path)\n"
    "          Input junctions tab-separated format.\n"
    "    -o, --output (std::filesystem::path)\n"
    "          The path of the vcf output file. If no path is given, will output to\n"
    "          standard output. Default: \"\".\n"
    "\n"
    "VERSION\n"
    "    Last update:\n"
    "    partitionJunctions version: 0.0.1\n"
    "    SeqAn version: 3.0.3\n"
    "\n"
    "LEGAL\n"
    "    Author: David Heller\n"
    "    SeqAn Copyright: 2006-2021 Knut Reinert, FU-Berlin; released under the\n"
    "    3-clause BSDL.\n"
};

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

TEST_F(find_deletions, help_page_argument)
{
    cli_test_result result = execute_app("find_deletions", "-h");
    std::string expected = help_page;

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(find_deletions, advanced_help_page_argument)
{
    cli_test_result result = execute_app("find_deletions", "-hh");
    std::string expected = help_page;

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
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
