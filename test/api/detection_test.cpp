#include <gtest/gtest.h>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>

#include "modules/sv_detection_methods/analyze_cigar_method.hpp"    // for the split read method
#include "modules/sv_detection_methods/analyze_sa_tag_method.hpp"   // for the cigar string method

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;

/* -------- detection methods tests -------- */

// TODO (irallia): implement test cases

TEST(junction_detection, cigar_string_simple_del)
{
    std::string const read_name = "read021";
    std::string const chromosome = "chr1";
    int32_t const query_start_pos = 1;
    std::vector<seqan3::cigar> cigar_string = {{5, 'S'_cigar_operation},
                                               {15, 'M'_cigar_operation},
                                               {6, 'D'_cigar_operation},
                                               {30, 'M'_cigar_operation}}; //5S15M6D30M
    seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5};

    // Deletion smaller than minimum variant size
    {
        std::vector<Junction> junctions_res{};
        uint64_t const min_var_length = 10;
        analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

        EXPECT_EQ(junctions_res.size(), 0);
    }

    // Deletion larger than minimum variant size
    {
        std::vector<Junction> junctions_res{};
        uint64_t const min_var_length = 5;
        analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

        Breakend new_breakend_1 {chromosome, 15, strand::forward};
        Breakend new_breakend_2 {chromosome, 22, strand::forward};
        std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, ""_dna5, read_name}};

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
        }
    }
}

TEST(junction_detection, cigar_string_del_padding)
{
    std::string const read_name = "read021";
    std::string const chromosome = "chr1";
    int32_t const query_start_pos = 1;
    std::vector<seqan3::cigar> cigar_string = {{5, 'S'_cigar_operation},
                                               {15, 'M'_cigar_operation},
                                               {60, 'P'_cigar_operation},
                                               {6, 'D'_cigar_operation},
                                               {30, 'M'_cigar_operation}}; //5S15M60P6D30M
    seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5};

    std::vector<Junction> junctions_res{};
    uint64_t const min_var_length = 5;
    analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

    Breakend new_breakend_1 {chromosome, 15, strand::forward};
    Breakend new_breakend_2 {chromosome, 22, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, ""_dna5, read_name}};

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
    }
}

TEST(junction_detection, cigar_string_simple_ins)
{
    std::string const read_name = "read021";
    std::string const chromosome = "chr1";
    int32_t const query_start_pos = 1;
    std::vector<seqan3::cigar> cigar_string = {{5, 'S'_cigar_operation},
                                               {9, 'M'_cigar_operation},
                                               {6, 'I'_cigar_operation},
                                               {30, 'M'_cigar_operation}}; //5S9M6I30M
    seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5};

    std::vector<Junction> junctions_res{};
    uint64_t const min_var_length = 5;
    analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

    Breakend new_breakend_1 {chromosome, 9, strand::forward};
    Breakend new_breakend_2 {chromosome, 10, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, "ATTTCG"_dna5, read_name}};

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
    }
}

TEST(junction_detection, cigar_string_ins_hardclip)
{
    std::string const read_name = "read021";
    std::string const chromosome = "chr1";
    int32_t const query_start_pos = 1;
    std::vector<seqan3::cigar> cigar_string = {{5, 'H'_cigar_operation},
                                               {9, 'M'_cigar_operation},
                                               {6, 'I'_cigar_operation},
                                               {35, 'M'_cigar_operation}}; //5H9M6I35M
    seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5};

    std::vector<Junction> junctions_res{};
    uint64_t const min_var_length = 5;
    analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

    Breakend new_breakend_1 {chromosome, 9, strand::forward};
    Breakend new_breakend_2 {chromosome, 10, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, "GATCGA"_dna5, read_name}};

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name());
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]);
    }
}

TEST(junction_detection, split_string)
{
    std::vector<std::string> test_strings{"a;a;a;aa", "b,  b,b,  bbb", "c c cccc", ";d;;d;"};
    {
        std::vector<std::string> resulting_strings{};
        split_string(test_strings[0], resulting_strings, ';');      // split "a;a;a;aa"     with ';'
        std::vector<std::string> expected{"a", "a", "a", "aa"};
        EXPECT_EQ(expected, resulting_strings);
    }
    {
        std::vector<std::string> resulting_strings{};
        split_string(test_strings[1], resulting_strings, ',');      // split "b,  b,b,  bbb" with ','
        std::vector<std::string> expected{ "b", "  b", "b", "  bbb" };
        EXPECT_EQ(expected, resulting_strings);
    }
    {
        std::vector<std::string> resulting_strings1{};
        split_string(test_strings[2], resulting_strings1, ' ');     // split "c c cccc"     with ' '
        std::vector<std::string> expected{"c", "c", "cccc"};
        EXPECT_EQ(expected, resulting_strings1);

        std::vector<std::string> resulting_strings2{};
        split_string(test_strings[2], resulting_strings2);          // split "c c cccc"     with default ' '
        EXPECT_EQ(resulting_strings1, resulting_strings2);
    }
    {
        std::vector<std::string> resulting_strings{};
        split_string(test_strings[3], resulting_strings, ';');      // split ";d;;d;" with ';'
        std::vector<std::string> expected{ "", "d", "", "d" };
        EXPECT_EQ(expected, resulting_strings);
    }
}

TEST(junction_detection, retrieve_aligned_segments)
{
    std::string const sa_tag = "chr1,100,+,6M94S,60,0;chr2,100,+,6S10M84S,60,0;chr1,106,+,16S10M74S,60,0;"
                               "chr1,116,-,60S14M26S,60,0;chr1,130,+,40S4M56S,60,0;chr1,150,+,44S6M50S,60,0;"
                               "chr1,156,+,90S10M,60,0;";
    std::vector<AlignedSegment> segments_res{};
    retrieve_aligned_segments(sa_tag, segments_res);

    AlignedSegment aligned_segment1 {strand::forward, "chr1", 100, 60, std::vector<seqan3::cigar>{{6, 'M'_cigar_operation},
                                                                                                  {94, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment2 {strand::forward, "chr2", 100, 60, std::vector<seqan3::cigar>{{6, 'S'_cigar_operation},
                                                                                                  {10, 'M'_cigar_operation},
                                                                                                  {84, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment3 {strand::forward, "chr1", 106, 60, std::vector<seqan3::cigar>{{16, 'S'_cigar_operation},
                                                                                                  {10, 'M'_cigar_operation},
                                                                                                  {74, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment4 {strand::reverse, "chr1", 116, 60, std::vector<seqan3::cigar>{{60, 'S'_cigar_operation},
                                                                                                  {14, 'M'_cigar_operation},
                                                                                                  {26, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment5 {strand::forward, "chr1", 130, 60, std::vector<seqan3::cigar>{{40, 'S'_cigar_operation},
                                                                                                  {4, 'M'_cigar_operation},
                                                                                                  {56, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment6 {strand::forward, "chr1", 150, 60, std::vector<seqan3::cigar>{{44, 'S'_cigar_operation},
                                                                                                  {6, 'M'_cigar_operation},
                                                                                                  {50, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment7 {strand::forward, "chr1", 156, 60, std::vector<seqan3::cigar>{{90, 'S'_cigar_operation},
                                                                                                  {10, 'M'_cigar_operation}}};
    std::vector<AlignedSegment> segments_expected_res{aligned_segment1,
                                                      aligned_segment2,
                                                      aligned_segment3,
                                                      aligned_segment4,
                                                      aligned_segment5,
                                                      aligned_segment6,
                                                      aligned_segment7};

    ASSERT_EQ(segments_expected_res.size(), segments_res.size());

    for (size_t i = 0; i < segments_expected_res.size(); ++i)
    {
        EXPECT_TRUE(segments_expected_res[i] == segments_res[i]);
    }
}

TEST(junction_detection, analyze_aligned_segments)
{
    AlignedSegment aligned_segment1 {strand::forward, "chr1", 100, 60, std::vector<seqan3::cigar>{{6, 'M'_cigar_operation},
                                                                                                  {94, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment2 {strand::forward, "chr2", 100, 60, std::vector<seqan3::cigar>{{6, 'S'_cigar_operation},
                                                                                                  {10, 'M'_cigar_operation},
                                                                                                  {84, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment3 {strand::forward, "chr1", 106, 60, std::vector<seqan3::cigar>{{16, 'S'_cigar_operation},
                                                                                                  {10, 'M'_cigar_operation},
                                                                                                  {74, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment4 {strand::reverse, "chr1", 116, 60, std::vector<seqan3::cigar>{{60, 'S'_cigar_operation},
                                                                                                  {14, 'M'_cigar_operation},
                                                                                                  {26, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment5 {strand::forward, "chr1", 130, 60, std::vector<seqan3::cigar>{{40, 'S'_cigar_operation},
                                                                                                  {4, 'M'_cigar_operation},
                                                                                                  {56, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment6 {strand::forward, "chr1", 150, 60, std::vector<seqan3::cigar>{{44, 'S'_cigar_operation},
                                                                                                  {6, 'M'_cigar_operation},
                                                                                                  {50, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment7 {strand::forward, "chr1", 156, 60, std::vector<seqan3::cigar>{{90, 'S'_cigar_operation},
                                                                                                  {10, 'M'_cigar_operation}}};
    std::vector<AlignedSegment> aligned_segments{aligned_segment1,
                                                 aligned_segment2,
                                                 aligned_segment3,
                                                 aligned_segment4,
                                                 aligned_segment5,
                                                 aligned_segment6,
                                                 aligned_segment7};

    // Minimum variant length = 10
    {
        std::vector<Junction> junctions_res{};
        seqan3::dna5_vector query_sequence = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"
                                              "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGACTGACAGGATA"_dna5};
        std::string const read_name = "read021";
        analyze_aligned_segments(aligned_segments,
                                 junctions_res,
                                 query_sequence,
                                 read_name,
                                 10,
                                 0);
        
        Breakend new_breakend_1 {"chr1", 105, strand::forward};
        Breakend new_breakend_2 {"chr2", 100, strand::forward};
        Breakend new_breakend_3 {"chr1", 106, strand::reverse};
        Breakend new_breakend_4 {"chr2", 109, strand::reverse};
        Breakend new_breakend_5 {"chr1", 115, strand::forward};
        Breakend new_breakend_6 {"chr1", 129, strand::reverse};
        Breakend new_breakend_7 {"chr1", 116, strand::reverse};
        Breakend new_breakend_8 {"chr1", 130, strand::forward};
        Breakend new_breakend_9 {"chr1", 133, strand::forward};
        Breakend new_breakend_10 {"chr1", 150, strand::forward};
        Breakend new_breakend_11 {"chr1", 155, strand::forward};
        Breakend new_breakend_12 {"chr1", 156, strand::forward};
        std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, ""_dna5, read_name},          //translocation
                                                     Junction{new_breakend_3, new_breakend_4, ""_dna5, read_name},          //translocation
                                                     Junction{new_breakend_5, new_breakend_6, ""_dna5, read_name},          //inversion
                                                     Junction{new_breakend_7, new_breakend_8, ""_dna5, read_name},          //inversion
                                                     Junction{new_breakend_9, new_breakend_10, ""_dna5, read_name},         //deletion
                                                     Junction{new_breakend_11,
                                                              new_breakend_12,
                                                              "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGAC"_dna5, read_name}}; //insertion

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]) << "Junction " << i << " unequal\nMate 1 equal: "
                                                                       << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                                                                       << "\nMate 2 equal: "
                                                                       << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                                                                       << "\n";
        }
    }

    // Minimum variant length = 20
    {
        std::vector<Junction> junctions_res{};
        seqan3::dna5_vector query_sequence = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"
                                              "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGACTGACAGGATA"_dna5};
        std::string const read_name = "read021";
        analyze_aligned_segments(aligned_segments,
                                 junctions_res,
                                 query_sequence,
                                 read_name,
                                 20,
                                 0);
        
        Breakend new_breakend_1 {"chr1", 105, strand::forward};
        Breakend new_breakend_2 {"chr2", 100, strand::forward};
        Breakend new_breakend_3 {"chr1", 106, strand::reverse};
        Breakend new_breakend_4 {"chr2", 109, strand::reverse};
        Breakend new_breakend_11 {"chr1", 155, strand::forward};
        Breakend new_breakend_12 {"chr1", 156, strand::forward};
        // The inversion and deletion are smaller than 20 bp and therefore not returned
        std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, ""_dna5, read_name},          //translocation
                                                     Junction{new_breakend_3, new_breakend_4, ""_dna5, read_name},          //translocation
                                                     Junction{new_breakend_11,
                                                              new_breakend_12,
                                                              "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGAC"_dna5, read_name}}; //insertion

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]) << "Junction " << i << " unequal\nMate 1 equal: "
                                                                       << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                                                                       << "\nMate 2 equal: "
                                                                       << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                                                                       << "\n";
        }
    }
}

TEST(junction_detection, analyze_sa_tag)
{
    // Primary alignment: chr1,116,-,10S14M26S,60,0;
    std::string const read_name = "read021";
    seqan3::sam_flag const flag{16u};
    std::string const chromosome = "chr1";
    int32_t const pos = 116;
    uint8_t const mapq = 60;
    std::vector<seqan3::cigar> test_cigar = {{10, 'S'_cigar_operation}, {14, 'M'_cigar_operation}, {26, 'S'_cigar_operation}};
    seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"
                               "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGACTGACAGGATA"_dna5};
    // Supplementary alignments
    std::string const sa_tag = "chr1,100,+,6M94S,60,0;chr2,100,+,6S10M84S,60,0;chr1,106,+,16S10M74S,60,0;"
                               "chr1,116,-,60S14M26S,60,0;chr1,130,+,40S4M56S,60,0;chr1,150,+,44S6M50S,60,0;"
                               "chr1,156,+,90S10M,60,0;";
    std::vector<Junction> junctions_res{};
    // Args
    cmd_arguments args{std::filesystem::path{},
                       std::filesystem::path{},
                       std::filesystem::path{},
                       std::vector<detection_methods>{cigar_string, split_read, read_pairs, read_depth},
                       simple_clustering,
                       sVirl_refinement_method,
                       10,
                       0};

    analyze_sa_tag(read_name, flag, chromosome, pos, mapq, test_cigar, seq, sa_tag, args, junctions_res);

    Breakend new_breakend_1 {"chr1", 105, strand::forward};
    Breakend new_breakend_2 {"chr2", 100, strand::forward};
    Breakend new_breakend_3 {"chr1", 106, strand::reverse};
    Breakend new_breakend_4 {"chr2", 109, strand::reverse};
    Breakend new_breakend_5 {"chr1", 115, strand::forward};
    Breakend new_breakend_6 {"chr1", 129, strand::reverse};
    Breakend new_breakend_7 {"chr1", 116, strand::reverse};
    Breakend new_breakend_8 {"chr1", 130, strand::forward};
    Breakend new_breakend_9 {"chr1", 133, strand::forward};
    Breakend new_breakend_10 {"chr1", 150, strand::forward};
    Breakend new_breakend_11 {"chr1", 155, strand::forward};
    Breakend new_breakend_12 {"chr1", 156, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2, ""_dna5, read_name},          //translocation
                                                 Junction{new_breakend_3, new_breakend_4, ""_dna5, read_name},          //translocation
                                                 Junction{new_breakend_5, new_breakend_6, ""_dna5, read_name},          //inversion
                                                 Junction{new_breakend_7, new_breakend_8, ""_dna5, read_name},          //inversion
                                                 Junction{new_breakend_9, new_breakend_10, ""_dna5, read_name},         //deletion
                                                 Junction{new_breakend_11,
                                                          new_breakend_12,
                                                          "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGAC"_dna5, read_name}}; //insertion

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(), junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]) << "Junction " << i << " unequal\nMate 1 equal: "
                                                                   << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                                                                   << "\nMate 2 equal: "
                                                                   << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                                                                   << "\n";
    }
}

// TODO (irallia): these tests should be implemented when the associated detection methods are implemented

// TEST(junction_detection, read_pairs_method_simple)
// {
//     testing::internal::CaptureStdout();

//     std::vector<Junction> resulting_junctions{};
//     analyze_read_pairs();

//     std::vector<Junction> expected_junctions{};
//     EXPECT_EQ(expected_junctions, resulting_junctions);
// }

// TEST(junction_detection, read_depth_method_simple)
// {
//     testing::internal::CaptureStdout();

//     std::vector<Junction> resulting_junctions{};
//     analyze_read_depth();

//     std::vector<Junction> expected_junctions{};
//     EXPECT_EQ(expected_junctions, resulting_junctions);
// }
