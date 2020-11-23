#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include <fstream>
#include <sstream>

#include "find_deletions/deletion_finding_and_printing.hpp"

TEST(group1, test_std_out)
{   
    find_and_print_deletions(DATADIR"detect_breakends_shorted.vcf", DATADIR"find_deletions_file_out_test.vcf");
    std::ifstream f;
    f.open(DATADIR"find_deletions_file_out_test.vcf");
    std::stringstream buffer;
    buffer << f.rdbuf();
    std::string expected
    {
        "##fileformat=VCFv4.2\n"
        "##source=iGenVarCaller\n"
        "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
        "chr21\t9435236\t.\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-50;END=9435286\n"
        "chr21\t11104574\t.\tN\t<DEL>\t60\tPASS\tSVTYPE=DEL;SVLEN=-248;END=11104822\n"
    };
    EXPECT_TRUE(f.is_open());
    EXPECT_EQ(buffer.str(), expected);
}
