#include <gtest/gtest.h>

#include <fstream>
#include <sstream>

#include "find_deletions/deletion_finding_and_printing.hpp"

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

TEST(group1, test_out_file)
{
    find_and_print_deletions(DATADIR"detect_breakends_shorted.vcf", DATADIR"find_deletions_file_out_test.vcf");
    std::ifstream f;
    f.open(DATADIR"find_deletions_file_out_test.vcf");
    std::stringstream buffer;
    buffer << f.rdbuf();
    EXPECT_TRUE(f.is_open());
    EXPECT_EQ(buffer.str(), expected);
}

TEST(group1, test_std_out)
{
    testing::internal::CaptureStdout();
    find_and_print_deletions(DATADIR"detect_breakends_shorted.vcf", "");
    std::string std_cout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(std_cout, expected);
}
