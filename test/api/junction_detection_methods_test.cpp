#include <gtest/gtest.h>

#include <seqan3/core/debug_stream.hpp>

#include "variant_detection/variant_detection.hpp"

using seqan3::operator""_dna5;

std::string const empty_alignment_short_reads_file_path = "";
std::string const default_alignment_long_reads_file_path = DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam";
std::filesystem::path empty_path{};
std::vector<detection_methods> const all_methods{cigar_string, split_read, read_pairs, read_depth};
constexpr uint64_t sv_default_length = 30;

// Explanation for the strings:
// chr21\t41972616\tForward\tchr21\t41972617\tForward\t1\t1681
// INS - Chromosome: chr21; Position: 41972615; Orientation: Forward
//       Chromosome: chr21; Position: 41972616; Orientation: Forward
//       Supporting Reads: 1
//       Average insertion size between two positions: 1681bp
// std::string expected_res_cigar
// {
//     "chr21\t41972615\tForward\tchr21\t41972616\tForward\t1\t1681\n"
// };

// chr22\t17458417\tForward\tchr21\t41972615\tForward\t1\t2
// BND - Chromosome: chr22; Position: 17458417; Orientation: Forward
//       Chromosome: chr21; Position: 41972615; Orientation: Forward
//       Supporting Reads: 1
//       Average insertion size between two positions: 2bp
// std::string expected_res_split
// {
//     "chr22\t17458417\tForward\tchr21\t41972615\tForward\t1\t2\n"
//     "chr22\t17458418\tForward\tchr21\t41972616\tForward\t2\t0\n"
// };

// std::string expected_res_pair = "";
// std::string expected_res_depth = "";

std::string empty_res
{
    "##fileformat=VCFv4.3\n"
    "##source=iGenVarCaller\n"
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
};

void check_output(std::string expected_res)
{
    std::string std_cout = testing::internal::GetCapturedStdout();
    EXPECT_EQ(expected_res, std_cout);

    return;
}

/* -------- detection methods tests -------- */
/* -------- single methods -------- */

TEST(junction_detection, single_method_only)
{
    for (detection_methods method : all_methods)
    {
        testing::internal::CaptureStdout();
        seqan3::debug_stream << "-----------------------------------------------------------------------\n"
                             << "Test Method: " << method << '\n'
                             << "-----------------------------------------------------------------------\n";
        detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                          default_alignment_long_reads_file_path,
                                          {method},
                                          simple_clustering,
                                          no_refinement,
                                          sv_default_length,
                                          empty_path);

        //TODO (eldariont): Currently, this (CLI-like) test compares the stdout of the method with an expected string in VCF format.
        //We need to replace this with real API tests that check inidividual parts of the pipeline.
        if (method == 0)
        {
            // check_output_and_cleanup(expected_res_cigar);
            check_output(empty_res);
        }
        else if (method == 1)
        {
            // check_output_and_cleanup(expected_res_split);
            check_output(empty_res);
        }
        else if (method == 2)
        {
            // check_output_and_cleanup(expected_res_pair);
            check_output(empty_res);
        }
        else // (method == 3)
        {
            // check_output_and_cleanup(expected_res_depth);
            check_output(empty_res);
        }
    }
}

/* -------- method pairs -------- */

TEST(junction_detection, method_pairs)
{
    for (detection_methods method_i : all_methods)
    {
        for (detection_methods method_j : all_methods)
        {
            if (method_i >= method_j) continue; // only check in one order to safe time

            testing::internal::CaptureStdout();
            seqan3::debug_stream << "-----------------------------------------------------------------------\n"
                                    << "Test Methods: " << method_i << ", " << method_j << '\n'
                                    << "-----------------------------------------------------------------------\n";
            detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                              default_alignment_long_reads_file_path,
                                              {method_i, method_j},
                                              simple_clustering,
                                              no_refinement,
                                              sv_default_length,
                                              empty_path);

            //TODO (eldariont): Currently, this (CLI-like) test compares the stdout of the method with an expected string in VCF format.
            //We need to replace this with real API tests that check inidividual parts of the pipeline.
            if (method_i == 0 && method_j == 1)
            {
                // check_output_and_cleanup(expected_res_cigar + expected_res_split);
                check_output(empty_res);
            }
            else if (method_i == 0 && method_j == 2)
            {
                // check_output_and_cleanup(expected_res_cigar + expected_res_pair);
                check_output(empty_res);
            }
            else if (method_i == 0 && method_j == 3)
            {
                // check_output_and_cleanup(expected_res_cigar + expected_res_depth);
                check_output(empty_res);
            }
            else if (method_i == 1 && method_j == 2)
            {
                // check_output_and_cleanup(expected_res_split + expected_res_pair);
                check_output(empty_res);
            }
            else if (method_i == 1 && method_j == 3)
            {
                // check_output_and_cleanup(expected_res_split + expected_res_depth);
                check_output(empty_res);
            }
            else // (method_i == 2 && method_j == 3)
            {
                // check_output_and_cleanup(expected_res_pair + expected_res_depth);
                check_output(empty_res);
            }
        }
    }
}

/* -------- method triples -------- */

TEST(junction_detection, method_triples)
{
    for (detection_methods method_i : all_methods)
    {
        for (detection_methods method_j : all_methods)
        {
            if (method_i >= method_j) continue; // only check in one order to safe time
            for (detection_methods method_k : all_methods)
            {
                if (method_j >= method_k) continue; // only check in one order to safe time

                testing::internal::CaptureStdout();
                seqan3::debug_stream << "-----------------------------------------------------------------------\n"
                                    << "Test Methods: " << method_i << ", " << method_j << ", " << method_k
                                    << '\n'
                                    << "-----------------------------------------------------------------------\n";
                detect_variants_in_alignment_file(empty_alignment_short_reads_file_path,
                                                  default_alignment_long_reads_file_path,
                                                  {method_i, method_j, method_k},
                                                  simple_clustering,
                                                  no_refinement,
                                                  sv_default_length,
                                                  empty_path);

                //TODO (eldariont): Currently, this (CLI-like) test compares the stdout of the method with an expected string in VCF format.
                //We need to replace this with real API tests that check inidividual parts of the pipeline.
                if (method_i == 0 && method_j == 1 && method_k == 2)
                {
                    // check_output_and_cleanup(expected_res_cigar + expected_res_split + expected_res_pair);
                    check_output(empty_res);
                }
                else if (method_i == 0 && method_j == 1 && method_k == 3)
                {
                    // check_output_and_cleanup(expected_res_cigar + expected_res_split + expected_res_depth);
                    check_output(empty_res);
                }
                else if (method_i == 0 && method_j == 2 && method_k == 3)
                {
                    // check_output_and_cleanup(expected_res_cigar + expected_res_pair + expected_res_depth);
                    check_output(empty_res);
                }
                else if (method_i == 1 && method_j == 2 && method_k == 3)
                {
                    // check_output_and_cleanup(expected_res_split + expected_res_pair + expected_res_depth);
                    check_output(empty_res);
                }
            }
        }
    }
}

TEST(junction_detection, detect_junctions_in_long_reads_sam_file)
{
    std::vector<Junction> junctions_res{};

    detect_junctions_in_long_reads_sam_file(junctions_res,
                                            default_alignment_long_reads_file_path,
                                            all_methods,
                                            simple_clustering,
                                            sVirl_refinement_method,
                                            sv_default_length);

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

    EXPECT_EQ(junctions_expected_res.size(), junctions_res.size());

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
