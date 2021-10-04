#include <gtest/gtest.h>

#include <fstream>

#include <seqan3/io/exception.hpp>

#include "iGenVar.hpp"                              // for global variable gVerbose
#include "variant_detection/variant_detection.hpp"  // for detect_junctions_in_long_reads_sam_file()

using seqan3::operator""_dna5;

std::string const default_alignment_short_reads_file_path = DATADIR"paired_end_mini_example.sam";
std::string const default_alignment_long_reads_file_path = DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam";
std::filesystem::path const empty_path{};
std::string default_vcf_sample_name{"MYSAMPLE"};
constexpr int16_t default_threads = 1;
std::vector<detection_methods> const default_methods{cigar_string, split_read, read_pairs, read_depth};
constexpr int32_t default_min_length = 30;
constexpr int32_t default_max_var_length = 1000000;
constexpr int32_t default_max_tol_inserted_length = 50;
constexpr int32_t default_max_tol_deleted_length = 50;
constexpr int32_t default_max_overlap = 10;
constexpr int32_t default_min_qual = 1;
constexpr int32_t default_partition_max_distance = 1000;
constexpr double default_hierarchical_clustering_cutoff = 0.5;

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
    std::map<std::string, int32_t> references_lengths{};

    cmd_arguments args{default_alignment_short_reads_file_path,
                       "",
                       empty_path, // empty genome path,
                       empty_path, // empty output path,
                       default_vcf_sample_name,
                       empty_path, // empty junctions path,
                       empty_path, // empty clusters path,
                       default_threads,
                       default_methods,
                       simple_clustering,
                       sVirl_refinement_method,
                       default_min_length,
                       default_max_var_length,
                       default_max_tol_inserted_length,
                       default_max_tol_deleted_length,
                       default_max_overlap,
                       default_min_qual,
                       default_partition_max_distance,
                       default_hierarchical_clustering_cutoff};
    detect_junctions_in_short_reads_sam_file(junctions_res, references_lengths, args);

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
    gVerbose = true;

    std::vector<Junction> junctions_res{};
    std::map<std::string, int32_t> references_lengths{};

    cmd_arguments args{"",
                       default_alignment_long_reads_file_path,
                       empty_path, // empty genome path,
                       empty_path, // empty output path,
                       default_vcf_sample_name,
                       empty_path, // empty junctions path,
                       empty_path, // empty clusters path,
                       default_threads,
                       default_methods,
                       simple_clustering,
                       sVirl_refinement_method,
                       default_min_length,
                       default_max_var_length,
                       default_max_tol_inserted_length,
                       default_max_tol_deleted_length,
                       default_max_overlap,
                       default_min_qual,
                       default_partition_max_distance,
                       default_hierarchical_clustering_cutoff};
    detect_junctions_in_long_reads_sam_file(junctions_res, references_lengths, args);

    std::string const chromosome_1 = "chr21";
    std::string const chromosome_2 = "chr22";
    int32_t const pos_ref_1 = 41972615;
    int32_t const pos_ref_2 = 41972616;
    int32_t const pos_ref_3 = 17458415;
    int32_t const pos_ref_4 = 41972615;
    int32_t const pos_ref_5 = 17458416;
    int32_t const pos_ref_6 = 41972616;
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
    Breakend new_breakend_2 {chromosome_1, pos_ref_2, strand::forward};
    Breakend new_breakend_5 {chromosome_2, pos_ref_3, strand::forward};
    Breakend new_breakend_6 {chromosome_1, pos_ref_4, strand::forward};
    Breakend new_breakend_7 {chromosome_2, pos_ref_5, strand::forward};
    Breakend new_breakend_8 {chromosome_1, pos_ref_6, strand::forward};

    size_t tandem_dup_count = 0;

    std::string const read_name_1 = "m2257/8161/CCS";
    std::string const read_name_2 = "m41327/11677/CCS";
    std::string const read_name_3 = "m21263/13017/CCS";
    std::string const read_name_4 = "m38637/7161/CCS";

    std::vector<Junction> junctions_expected_res
    {
        Junction{new_breakend_1, new_breakend_2, insertion_sequence_1, tandem_dup_count, read_name_1},
        Junction{new_breakend_5, new_breakend_6, "TA"_dna5, tandem_dup_count, read_name_2},
        Junction{new_breakend_7, new_breakend_8, ""_dna5, tandem_dup_count, read_name_3},
        Junction{new_breakend_7, new_breakend_8, ""_dna5, tandem_dup_count, read_name_4}
    };

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
        // For debugging #include <seqan3/core/debug_stream.hpp> and use:
        // seqan3::debug_stream << "-----------------------------------------------------------------------------------\n"
        //                      << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1()) << ":  "
        //                      << junctions_expected_res[i].get_mate1() << "\n == " << junctions_res[i].get_mate1() << "\n"
        //                      << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2()) << ":  "
        //                      << junctions_expected_res[i].get_mate2() << "\n == " << junctions_res[i].get_mate2() << "\n"
        //                      << junctions_expected_res[i].get_inserted_sequence() << " == " << junctions_res[i].get_inserted_sequence() << "\n"
        //                      << junctions_expected_res[i].get_tandem_dup_count() << " == " << junctions_res[i].get_tandem_dup_count() << "\n";
    }
}

TEST(input_file, long_read_sam_file_unsorted)
{
    std::vector<Junction> junctions_res{};
    std::map<std::string, int32_t> references_lengths{};

    // Create a blank SAM file without a sorting indicator.
    std::filesystem::path const tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory
    std::filesystem::path unsorted_sam_path{tmp_dir/"unsorted.sam"};
    std::ofstream unsorted_sam{unsorted_sam_path.c_str()};
    unsorted_sam << "@HD\tVN:1.6\n"
                 << "@SQ\tSN:testchr\tLN:1000\n"
                 << "test1\t16\ttestchr\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n";
    unsorted_sam.close();

    cmd_arguments args{"",
                       unsorted_sam_path,
                       empty_path, // empty genome path,
                       empty_path, // empty output path,
                       default_vcf_sample_name,
                       empty_path, // empty junctions path,
                       empty_path, // empty clusters path,
                       default_threads,
                       default_methods,
                       simple_clustering,
                       no_refinement,
                       default_min_length,
                       default_max_var_length,
                       default_max_tol_inserted_length,
                       default_max_tol_deleted_length,
                       default_max_overlap,
                       default_min_qual,
                       default_partition_max_distance,
                       default_hierarchical_clustering_cutoff};
    EXPECT_THROW(detect_junctions_in_long_reads_sam_file(junctions_res,
                                                         references_lengths,
                                                         args), seqan3::format_error);

    std::filesystem::remove(unsorted_sam_path);
}

TEST(input_file, short_and_long_read_sam_file_with_different_references_lengths)
{
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    std::string const expected_err
    {
        "The cigar string method for short reads is not yet implemented.\n"
        "The split read method for short reads is not yet implemented.\n"
        "The read depth method for short reads is not yet implemented.\n"
        "The cigar string method for short reads is not yet implemented.\n"
        "The split read method for short reads is not yet implemented.\n"
        "The read depth method for short reads is not yet implemented.\n"
        "The cigar string method for short reads is not yet implemented.\n"
        "The split read method for short reads is not yet implemented.\n"
        "The read depth method for short reads is not yet implemented.\n"
        "Warning: The reference id chr2 was found twice in the input files with different length: 1001 and 1005\n"
        "Warning: The reference id chr4 was found twice in the input files with different length: 1004 and 1005\n"
        "The read depth method for long reads is not yet implemented.\n"
        "The read depth method for long reads is not yet implemented.\n"
        "The read depth method for long reads is not yet implemented.\n"
    };

    std::vector<Junction> junctions_res{};
    std::map<std::string, int32_t> references_lengths{};

    std::filesystem::path const tmp_dir = std::filesystem::temp_directory_path();     // get the temp directory

    // Create a blank short read SAM file with SQ header tag with different length of one reference.
    std::filesystem::path short_sam_path{tmp_dir/"short.sam"};
    std::ofstream short_sam{short_sam_path.c_str()};
    short_sam << "@HD\tVN:1.6\tSO:coordinate\n"
              << "@SQ\tSN:chr1\tLN:1000\n"          // chr1 present in both files with same length
              << "@SQ\tSN:chr2\tLN:1001\n"          // chr2 present in both files with different length -> Warning
              << "@SQ\tSN:chr3\tLN:1002\n"          // chr3 present in same file twice with same length
              << "@SQ\tSN:chr3\tLN:1002\n"
              << "test1\t16\tchr1\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n"
              << "test1\t16\tchr2\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n"
              << "test1\t16\tchr3\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n";
    short_sam.close();

    // Create a blank long read SAM file with SQ header tag with different length of one reference.
    std::filesystem::path long_sam_path{tmp_dir/"long.sam"};
    std::ofstream long_sam{long_sam_path.c_str()};
    long_sam << "@HD\tVN:1.6\tSO:coordinate\n"
             << "@SQ\tSN:chr1\tLN:1000\n"
             << "@SQ\tSN:chr2\tLN:1005\n"
             << "@SQ\tSN:chr4\tLN:1004\n"           // chr4 present in same file twice with different length -> Warning
             << "@SQ\tSN:chr4\tLN:1005\n"
             << "test1\t16\tchr1\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n"
             << "test1\t16\tchr2\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n"
             << "test1\t16\tchr4\t1\t60\t10M\t=\t1\t0\tGCGCGCGCGC\tFFFFFFFFFF\n";
    long_sam.close();

    cmd_arguments args{short_sam_path,
                       long_sam_path,
                       empty_path, // empty genome path
                       empty_path, // empty output path
                       default_vcf_sample_name,
                       empty_path, // empty junctions path
                       empty_path, // empty clusters path
                       default_threads,
                       default_methods,
                       simple_clustering,
                       no_refinement,
                       default_min_length,
                       default_max_var_length,
                       default_max_tol_inserted_length,
                       default_max_tol_deleted_length,
                       default_max_overlap,
                       default_min_qual,
                       default_partition_max_distance,
                       default_hierarchical_clustering_cutoff};
    EXPECT_NO_THROW(detect_junctions_in_short_reads_sam_file(junctions_res, references_lengths, args));
    EXPECT_NO_THROW(detect_junctions_in_long_reads_sam_file(junctions_res, references_lengths, args));

    std::string result_out = testing::internal::GetCapturedStdout();
    EXPECT_EQ("", result_out);
    std::string result_err = testing::internal::GetCapturedStderr();
    EXPECT_EQ(expected_err, result_err);

    std::filesystem::remove(short_sam_path);
    std::filesystem::remove(long_sam_path);
}
