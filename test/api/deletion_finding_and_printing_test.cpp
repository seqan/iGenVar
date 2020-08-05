#include <gtest/gtest.h>
#include <gmock/gmock.h>

#include "find_deletions/deletion_finding_and_printing.hpp"

TEST(group1, test_std_out)
{
    std::string expected{"##fileformat=VCFv4.2\n"
                         "##source=iGenVarCaller\n"
                         "CHROM	POS	ID	REF	ALT	QUAL	FILTER	INFO\n"
                         "chr21	9435236	.	N	<DEL>	60	PASS	SVTYPE=DEL;SVLEN=-50;END=9435286\n"
                         "chr21	11104574	.	N	<DEL>	60	PASS	SVTYPE=DEL;SVLEN=-248;END=11104822\n"};

    testing::internal::CaptureStdout();
    find_and_print_deletions(DATADIR"detect_breakends_shorted.vcf");
    std::string std_cout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(expected, std_cout);
}
