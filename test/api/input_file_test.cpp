#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/io/exception.hpp>

#include "variant_detection/variant_detection.hpp"  // for detect_junctions_in_long_reads_sam_file()

using seqan3::operator""_dna5;

std::string const default_alignment_short_reads_file_path = DATADIR"paired_end_mini_example.sam";
std::string const default_alignment_long_reads_file_path = DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam";
std::filesystem::path const empty_output_path{};
std::vector<detection_methods> const default_methods{cigar_string, split_read, read_pairs, read_depth};
constexpr uint64_t sv_default_length = 30;

// Explanation for the strings:
// chr21\t41972615\tForward\tchr21\t41972616\tForward\t1\t1681
// INS - Chromosome: chr21; Position: 41972615; Orientation: Forward
//       Chromosome: chr21; Position: 41972616; Orientation: Forward
//       Supporting Reads: 1
//       Average insertion size between two positions: 1681bp
// chr22\t17458417\tForward\tchr21\t41972615\tForward\t1\t2
// BND - Chromosome: chr22; Position: 17458417; Orientation: Forward
//       Chromosome: chr21; Position: 41972615; Orientation: Forward
//       Supporting Reads: 1
//       Average insertion size between two positions: 2bp
// std::string expected_res
// {
//     "chr21\t41972615\tForward\tchr21\t41972616\tForward\t1\t1681\n"
//     "chr22\t17458417\tForward\tchr21\t41972615\tForward\t1\t2\n"
//     "chr22\t17458418\tForward\tchr21\t41972616\tForward\t2\t0\n"
// };

std::string empty_res
{
    "##fileformat=VCFv4.3\n"
    "##source=iGenVarCaller\n"
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
};

TEST(input_file, detect_junctions_in_short_read_sam_file)
{
    std::vector<Junction> junctions_res{};

    cmd_arguments args{default_alignment_short_reads_file_path,
                       "",
                       empty_output_path,
                       default_methods,
                       simple_clustering,
                       sVirl_refinement_method,
                       sv_default_length};
    detect_junctions_in_short_reads_sam_file(junctions_res,
                                             args);

    std::vector<Junction> junctions_expected_res{};

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
        // For debugging #include <seqan3/core/debug_stream.hpp> and use:
        // seqan3::debug_stream << "-----------------------------------------------------------------------------------\n"
        //                      << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1()) << ": \n"
        //                      << junctions_expected_res[i].get_mate1() << " == " << junctions_res[i].get_mate1() << "\n"
        //                      << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2()) << ": \n"
        //                      << junctions_expected_res[i].get_mate2() << " == " << junctions_res[i].get_mate2() << "\n";
    }
}

TEST(input_file, detect_junctions_in_long_reads_sam_file)
{
    std::vector<Junction> junctions_res{};

    cmd_arguments args{"",
                       default_alignment_long_reads_file_path,
                       empty_output_path,
                       default_methods,
                       simple_clustering,
                       sVirl_refinement_method,
                       sv_default_length};
    detect_junctions_in_long_reads_sam_file(junctions_res,
                                            args);

    std::string const chromosome_1 = "chr21";
    std::string const chromosome_2 = "chr22";
    int32_t const pos_ref_1 = 41972615;
    int32_t const pos_ref_2 = 17458417;
    int32_t const pos_ref_3 = 41972615;
    int32_t const pos_ref_4 = 17458418;
    int32_t const pos_ref_5 = 41972616;
    // std::string const insertion_allele_id_1 = "0";

    seqan3::dna5_vector const insertion_sequence_1 = "GAGTGGACCTCAGCAAACTCCCAGTAGAGCTGCAGCAGAGGGGGTCTGAC"
                                                     "TGTTAGGATGGAAAACTAACAAACAGAAGGCAATAGCATCAACAACAACA"
                                                     "AAAAAAAACTGCCACACAAAAACCGCTATCCGAAGATCACCAACATCAAA"
                                                     "GATCGAACAGGTAGACAAATGACGAAGAGAGGAAAAAACAGTGCAAAAAA"
                                                     "GGCTGAAAAGCCCAGAACCACCTCTTCCCTCCAGAGGATCACAACTCCTC"
                                                     "ACCAGGAAAGGGAACAAAACTGCACAGAGAAGAGTTTGACCAATGACAGA"
                                                     "AGTAGGCTTCAGCAGAATGGGTAATAACTCCTCTGAGCTAAAGGAGCATG"
                                                     "TTCCTACCCTAATGCAAGGAAGCTAAAGATACTTGATAAAAGGTTACAGG"
                                                     "AACTGCTAACTAAATAAACCAGTTCAGAGAAGAACATAAATGACCTAAAT"
                                                     "GGAGCTGAAAAACACAGCATGAGAACTTCATGAAGGAATACACAAGGTAT"
                                                     "CAACAGACAAGTCCATCAGGCAGAAGAAAGGGATATCAGAGATTGAAGAT"
                                                     "CAACTTAATGAAATAAAGCATGCAAGACGAGATTAGTGAGAAAAAAGAAT"
                                                     "TAAAAGAAATGAGCAAAGGCCTCAAGGAAATATGGGACTATGTGTAAAAG"
                                                     "ACCAAGCATACGGTTTGATTGGTGTATGTGAAAATGACAGGGAAAATGGA"
                                                     "ACCAAGTTGGAAAACACTCTTCAGGATATCATGCAGGAGAACCTCCCAAC"
                                                     "CTAGCAAGAGAAGCCAACATTCACATTCAGGAAATACAAGAGAACACCAC"
                                                     "CAAGATACTCCTTGAGAAGTAGCAAACCCCCAAGACACATAATTGTTCAG"
                                                     "ATTCAGGCAAGGGTGAAAATGAAGGAAAAAATGCTAAGAGAGCCAGAGAG"
                                                     "AAAGGTATGGGTTATCCACAAAAGGGCCAGCCATCAGACTAAGAGCAATT"
                                                     "CTCTGCAGAAACCCTACAACCAGAAGAGAGAAGGGGCCAATATTCAACAT"
                                                     "TCTTAAAGAAAAGAATTTTCAACCCAGAATTTCATATCCAGCCAAAACAA"
                                                     "AGCTTCGTAAGTGAAGGAGAAATAAATTTCTTTACAGACAAGCAAATGAC"
                                                     "TGAAGAAGATTTTTGTCACCACCATGCCTGCCTTACAAGATCTCCTGAAG"
                                                     "GAAGCACTAAGACATGGGAAGGAAAAAATCCAGTACCAAGCCACTGCTAA"
                                                     "ACCATACCAAAATGTAGAGACTCAATGCTTAGGATAGGAAACTGCATCAA"
                                                     "CTAGCAGTCAAAATAACCAGCTAGCATTCATAATGACAGGATCAAATTCA"
                                                     "GACCACATACAATTATTAACCTTAAATGTAAATGGGCTAAATGCCGCAAT"
                                                     "TAAAAGACACATCACTGGCAAATTGGATAAAGAGTCAAGCCCAATCGGTG"
                                                     "TGCTGTTATTCAAGGAGACACCACTCTCACGTGCAAGAGACACAGATAGG"
                                                     "CTCGAAAATGAATAGGGATGAAGGAAGATTACCAAGCAAATGGAAAGCAA"
                                                     "AAAAAAAAAAAGCAGGGGTTGCAAATCCTAGTCTCTGTAAAACACTACTT"
                                                     "TAAACCAAGAAAGATCAAAAGAGACAAAGAGGGCTTTAATAATGGTAATA"
                                                     "GGGGGATTAATTCAACAAGAAGAGTTAACTATCCTAAATATATATGCTGC"
                                                     "CTAATACAGGCACACCCAGATTCATAAAGCA"_dna5;

    Breakend new_breakend_1 {chromosome_1, pos_ref_1, strand::forward};
    Breakend new_breakend_2 {chromosome_1, pos_ref_1+1, strand::forward};
    Breakend new_breakend_5 {chromosome_2, pos_ref_2, strand::forward};
    Breakend new_breakend_6 {chromosome_1, pos_ref_3, strand::forward};
    Breakend new_breakend_7 {chromosome_2, pos_ref_4, strand::forward};
    Breakend new_breakend_8 {chromosome_1, pos_ref_5, strand::forward};
    Breakend new_breakend_9 {chromosome_1, pos_ref_5, strand::forward};

    std::string const read_name_1 = "m2257/8161/CCS";
    std::string const read_name_2 = "m41327/11677/CCS";
    std::string const read_name_3 = "m21263/13017/CCS";
    std::string const read_name_4 = "m38637/7161/CCS";

    std::vector<Junction> junctions_expected_res
    {   Junction{new_breakend_1, new_breakend_2, insertion_sequence_1, read_name_1},
        Junction{new_breakend_5, new_breakend_6, "TA"_dna5, read_name_2},
        Junction{new_breakend_7, new_breakend_8, ""_dna5, read_name_3},
        Junction{new_breakend_7, new_breakend_9, ""_dna5, read_name_4}
    };

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
        // For debugging #include <seqan3/core/debug_stream.hpp> and use:
        // seqan3::debug_stream << "-----------------------------------------------------------------------------------\n"
        //                      << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1()) << ": \n"
        //                      << junctions_expected_res[i].get_mate1() << " == " << junctions_res[i].get_mate1() << "\n"
        //                      << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2()) << ": \n"
        //                      << junctions_expected_res[i].get_mate2() << " == " << junctions_res[i].get_mate2() << "\n";
    }
}

TEST(input_file, long_read_sam_file_unsorted)
{
    std::vector<Junction> junctions_res{};
    // Create a blank SAM file without a sorting indicator.
    std::filesystem::path const tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::path unsorted_sam_path{tmp_dir/"unsorted.sam"};
    std::ofstream unsorted_sam{unsorted_sam_path.c_str()};
    unsorted_sam << "@HD\tVN:1.6\n" <<
                    "@SQ\tSN:testchr\tLN:1000\n" <<
                    "test1\t16\ttestchr\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n";
    unsorted_sam.close();

    cmd_arguments args{"",
                       unsorted_sam_path,
                       empty_output_path,
                       default_methods,
                       simple_clustering,
                       no_refinement,
                       sv_default_length};
    EXPECT_THROW(detect_junctions_in_long_reads_sam_file(junctions_res,
                                                         args), seqan3::format_error);

    std::filesystem::remove(unsorted_sam_path);
}
