#include "api_test.hpp"

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/io/sam_file/sam_flag.hpp>

#include "modules/sv_detection_methods/analyze_cigar_method.hpp"        // for the split read method
#include "modules/sv_detection_methods/analyze_split_read_method.hpp"   // for the cigar string method

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;

/* -------- detection methods tests -------- */

// TODO (irallia): implement test cases <- (23.7.21, irallia) which cases are done / still open?

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
        int32_t const min_var_length = 10;
        analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

        EXPECT_EQ(junctions_res.size(), 0u);
    }

    // Deletion larger than minimum variant size
    {
        std::vector<Junction> junctions_res{};
        int32_t const min_var_length = 5;
        analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

        Breakend new_breakend_1 {chromosome, 15, strand::forward};
        Breakend new_breakend_2 {chromosome, 22, strand::forward};
        std::vector<Junction> junctions_expected_res{Junction{new_breakend_1,
                                                              new_breakend_2,
                                                              ""_dna5,
                                                              0,
                                                              read_name}};

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
    int32_t const min_var_length = 5;
    analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

    Breakend new_breakend_1 {chromosome, 15, strand::forward};
    Breakend new_breakend_2 {chromosome, 22, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1,
                                                          new_breakend_2,
                                                          ""_dna5,
                                                          0,
                                                          read_name}};

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
    int32_t const min_var_length = 5;
    analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

    Breakend new_breakend_1 {chromosome, 9, strand::forward};
    Breakend new_breakend_2 {chromosome, 10, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1,
                                                          new_breakend_2,
                                                          "ATTTCG"_dna5,
                                                          0,
                                                          read_name}};

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
    int32_t const min_var_length = 5;
    analyze_cigar(read_name, chromosome, query_start_pos, cigar_string, seq, junctions_res, min_var_length);

    Breakend new_breakend_1 {chromosome, 9, strand::forward};
    Breakend new_breakend_2 {chromosome, 10, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1,
                                                          new_breakend_2,
                                                          "GATCGA"_dna5,
                                                          0,
                                                          read_name}};

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
    std::string const sa_tag = "chr1,101,+,6M94S,60,0;chr2,101,+,6S10M84S,60,0;chr1,107,+,16S10M74S,60,0;"
                               "chr1,117,-,60S14M26S,60,0;chr1,131,+,40S4M56S,60,0;chr1,151,+,44S6M50S,60,0;"
                               "chr1,157,+,90S10M,60,0;";
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
    auto verboseGuard = verbose_guard(true); // will reset back to the original state, after leaving this scope

    // read: |as_1||..as_2..||..as_3..||....as_4....||_5||as_6||..as_7..|
    //
    //                  100   106       116           130                 150   156
    //                   |     |         |             |                   |     |
    // chr1:             .....................................................................
    // aligned_segment1: |as_1|
    // aligned_segment3:       |..as_3..|
    // aligned_segment4:                 |....as_4....|
    // aligned_segment5:                               |_5|
    // aligned_segment6:                                                   |as_6|
    // aligned_segment7:                                                         |..as_7..|

    //                  100
    //                   |
    // chr2:             .....................................................................
    // aligned_segment2: |..as_2..|

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
        std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2,
                                                              ""_dna5,
                                                              0, read_name},         // translocation
                                                     Junction{new_breakend_3, new_breakend_4,
                                                              ""_dna5,
                                                              0, read_name},         // translocation
                                                     Junction{new_breakend_5, new_breakend_6,
                                                              ""_dna5,
                                                              0, read_name},         // inversion
                                                     Junction{new_breakend_7, new_breakend_8,
                                                              ""_dna5,
                                                              0, read_name},         // inversion
                                                     Junction{new_breakend_9, new_breakend_10,
                                                              ""_dna5,
                                                              0, read_name},         // deletion
                                                     Junction{new_breakend_11, new_breakend_12,
                                                              "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGAC"_dna5,
                                                              0, read_name}};        // insertion

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(),
                      junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
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
        std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2,
                                                              ""_dna5,
                                                              0, read_name},         // translocation
                                                     Junction{new_breakend_3, new_breakend_4,
                                                              ""_dna5,
                                                              0, read_name},         // translocation
                                                     Junction{new_breakend_11,
                                                              new_breakend_12,
                                                              "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGAC"_dna5,
                                                              0, read_name}};        // insertion

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(),
                      junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]) << "Junction " << i << " unequal\nMate 1 equal: "
                                                                       << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                                                                       << "\nMate 2 equal: "
                                                                       << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                                                                       << "\n";
        }
    }
}

TEST(junction_detection, overlapping_segments)
{
    AlignedSegment aligned_segment1 {strand::forward, "chr1", 100, 60, std::vector<seqan3::cigar>{{20, 'M'_cigar_operation},
                                                                                                  {30, 'S'_cigar_operation}}};
    AlignedSegment aligned_segment2 {strand::forward, "chr1", 200, 60, std::vector<seqan3::cigar>{{15, 'S'_cigar_operation},
                                                                                                  {35, 'M'_cigar_operation}}};
    std::vector<AlignedSegment> aligned_segments{aligned_segment1, aligned_segment2};

    std::vector<Junction> junctions_res{};
    seqan3::dna5_vector query_sequence = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5};
    std::string const read_name = "read021";
    ASSERT_NO_THROW(analyze_aligned_segments(aligned_segments,
                                             junctions_res,
                                             query_sequence,
                                             read_name,
                                             10,
                                             10));

    // Deletion from two overlapping alignment segments (overlap of 5bp)
    Breakend new_breakend_1 {"chr1", 119, strand::forward};
    Breakend new_breakend_2 {"chr1", 205, strand::forward};
    std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2,
                                                          ""_dna5,
                                                          0, read_name}};

    ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

    for (size_t i = 0; i < junctions_expected_res.size(); ++i)
    {
        EXPECT_EQ(junctions_expected_res[i].get_read_name(),
                  junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
        EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i]) << "Junction " << i << " unequal\nMate 1 equal: "
                                                                   << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                                                                   << "\nMate 2 equal: "
                                                                   << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                                                                   << "\n";
    }
}

TEST(junction_detection, analyze_sa_tag)
{
    auto verboseGuard = verbose_guard(true); // will reset back to the original state, after leaving this scope

    // Args
    cmd_arguments args{std::filesystem::path{}, // alignment_short_reads_file_path
                       std::filesystem::path{}, // alignment_long_reads_file_path
                       std::filesystem::path{}, // genome_file_path
                       std::filesystem::path{}, // output_file_path
                       "MYSAMPLE",              // vcf_sample_name
                       std::filesystem::path{}, // junctions_file_path
                       std::filesystem::path{}, // clusters_file_path
                       1,                       // threads
                       std::vector<detection_methods>{cigar_string, split_read, read_pairs, read_depth},
                       simple_clustering,
                       sVirl_refinement_method,
                       9,                       // min_var_length,
                       1000000,                 // default max_var_length,
                       50,                      // default max_tol_inserted_length,
                       50,                      // default max_tol_deleted_length,
                       0,                       // max_overlap,
                       1,                       // default min_qual,
                       1000,                    // default partition_max_distance
                       0.5};                    // default hierarchical_clustering_cutoff

    { // Example 1
        seqan3::debug_stream << "----------------------------------First Example:----------------------------------\n";
        std::vector<Junction> junctions_res{};

        // pos  11        21  31       41  51                        61        71  81      91         101 111    121         131 141     151
        //      |         |   |        |   |                         |         |   |       |           |   |      |           |   |       |
        // chr1 -1->      -2->-------> -3->------------>             -4->      -5->------->-6->        -7->------>-8->        -9->------->-10->
        //      ||||      |||||||||||| |||||||||||||||||\\\\\\\\\\\\ ||||      ||||||||||||||||        ||||       ||||        ||||
        // read --->-TRA->---><--INV-- --->-DUP:TANDEM->-DUP:TANDEM->--->-TRA->--->-DUP_1->--->--INS-->--->..DEL..--->-DUP_2->--->
        //          ||||||                                               ||||||                                       ||||||||
        // chr2 --->----->--->                                  chr1 -9->----->-10->                         chr1 -5->------->-6->
        //          |     |          | |   |            |            |   |     |   |       |   |       |          |   |       |
        // read_pos 11    21        31 41  51           61           71  81    91 101     111 121     131        141 151     161

        // Primary alignment: chr1,71,+,90S30M49S,60,0;
        std::string read_name = "read";
        seqan3::sam_flag flag{0};
        std::string chromosome = "chr1";
        int32_t pos = 70; // Decrement by 1 because scanned position is already 0-based
        uint8_t mapq = 60;
        std::vector<seqan3::cigar> test_cigar = {{90, 'S'_cigar_operation},
                                                 {30, 'M'_cigar_operation},
                                                 {49, 'S'_cigar_operation}};
        seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC" //   1 -  50
                                   "GCGATACGCGTCGCAACTACGACGCCATCAGCAGGCGACTGACAGGATAT" //  51 - 100
                                   "GGGTTAGTCCCTATCGATCGTACGTACAGATTCGATGGGGGAATTTGGAT" // 101 - 150
                                   "ACGGTTACGGGAGACCCTGA"_dna5};                        // 151 - 170
        // Supplementary alignments
        std::string sa_tag = "chr1,11,+,10M159S,60,0;"
                             "chr2,102,+,10S10M149S,60,0;"  // TRA 1 (interspersed): chr1 20 -> chr2 102 & chr2 111 -> chr1 21
                             "chr1,21,+,20S10M139S,60,0;"
                             "chr1,31,-,129S10M30S,60,0;"   // INV: (30,40] deleted, [31,41) inserted
                             "chr1,41,+,40S20M109S,60,0;"   // DUP:TANDEM: [51,60]
                             "chr1,51,+,60S20M89S,60,0;"    // DUP:TANDEM: [51,60]
                             "chr1,141,+,80S10M79S,60,0;"   // TRA 2 (intrachromosomal): 70 -> 141 & 150 -> 71
        // AS=30 -> primary  "chr1,71,+,90S30M49S,60,0;"    // DUP_1: matches to [81,90]
                                                            // INS: inserted after 100
                             "chr1,101,+,130S10M29S,60,0;"
                                                            // DEL: (110,120] deleted
                             "chr1,121,+,140S10M19S,60,0;"
                             "chr1,81,+,150S10M9S,60,0;"    // DUP_2: matches to [81,90]
                             "chr1,131,+,160S9M,60,0;";

        analyze_sa_tag(read_name, flag, chromosome, pos, mapq, test_cigar, seq, sa_tag, args, junctions_res);

        std::vector<Junction> junctions_expected_res =
        {                                                                     // +1 as we are 0 based but SAM is 1 based
            Junction{Breakend{"chr1", 19, strand::forward}, Breakend{"chr2", 101, strand::forward}, // TRA 1 (interspersed)
                     ""_dna5, 0, read_name},                                                        // chr1 20 -> chr2 102
            Junction{Breakend{"chr1", 20, strand::reverse}, Breakend{"chr2", 110, strand::reverse}, // TRA 1 (interspersed)
                     ""_dna5, 0, read_name},                                                        // chr2 111 -> chr1 21
            Junction{Breakend{"chr1", 29, strand::forward}, Breakend{"chr1", 39, strand::reverse},  // INV
                     ""_dna5, 0, read_name},                                                        // (30,40] deleted
            Junction{Breakend{"chr1", 30, strand::reverse}, Breakend{"chr1", 40, strand::forward},  // INV
                     ""_dna5, 0, read_name},                                                        // [31,41) inserted
            Junction{Breakend{"chr1", 50, strand::forward}, Breakend{"chr1", 59, strand::forward},  // DUP:TANDEM
                     "TCGCAACTAC"_dna5, 2, read_name},                                              // [51,60] duplicated
            Junction{Breakend{"chr1", 69, strand::forward}, Breakend{"chr1", 140, strand::forward}, // TRA 2 (intrachr.)
                     ""_dna5, 0, read_name},                                                        // behind 70 -> 141
            Junction{Breakend{"chr1", 70, strand::reverse}, Breakend{"chr1", 149, strand::reverse}, // TRA 2 (intrachr.)
                     ""_dna5, 0, read_name},                                                        // behind 150 -> 71
            Junction{Breakend{"chr1", 99, strand::forward}, Breakend{"chr1", 100, strand::forward}, // INS
                     "TACGTACAGA"_dna5, 0, read_name},                                              // inserted after 100
            Junction{Breakend{"chr1", 109, strand::forward}, Breakend{"chr1", 120, strand::forward},// DEL
                     ""_dna5, 0, read_name},                                                        // (110,120] deleted
            Junction{Breakend{"chr1", 80, strand::reverse}, Breakend{"chr1", 129, strand::reverse}, // DUP_2
                     ""_dna5, 0, read_name},                                                        // matches to [81,90]?!
            Junction{Breakend{"chr1", 89, strand::forward}, Breakend{"chr1", 130, strand::forward}, // DUP_2
                     ""_dna5, 0, read_name}                                                         // matches to [81,90]?!
        };

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(),
                      junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i])
                        << "Junction " << i << " unequal\n"
                        << "Mate 1 equal: " << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                        << "\nMate 2 equal: " << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                        << "\n";
        }

        // For debugging use:
        // print_compare_junction_vectors(junctions_expected_res, junctions_res);
    }

    { // Example 2 - multiple tandem duplication on read
        seqan3::debug_stream << "----------------------------------Second Example:----------------------------------\n";
        std::vector<Junction> junctions_res{};

        // pos  11  21                                                  31
        //      |   |                                                   |
        // chr1 -1->-----2------>                                       -3->
        //      |||||||||||||||||\\\\\\\\\\\\ \\\\\\\\\\\\ \\\\\\\\\\\\ ||||
        // read --->-DUP:TANDEM->-DUP:TANDEM->-DUP:TANDEM->-DUP:TANDEM->--->
        //          |            |            |            |            |
        // read_pos 11           21           31           41           51

        // Primary alignment: chr1,11,+,20M39S,60,0;
        std::string read_name = "read";
        seqan3::sam_flag flag{0};
        std::string chromosome = "chr1";
        int32_t pos = 10; // Decrement by 1 because scanned position is already 0-based
        uint8_t mapq = 60;
        std::vector<seqan3::cigar> test_cigar = {{20, 'M'_cigar_operation}, {39, 'S'_cigar_operation}};
        seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5}; // 1 - 50

        // Supplementary alignments
        std::string sa_tag = /*"chr1,11,+,20M39S,60,0;"*/   // AS=30 -> primary;
                             "chr1,21,+,20S10M29S,60,0;"    // DUP:TANDEM: [21,30]
                             "chr1,21,+,30S10M19S,60,0;"    // DUP:TANDEM: [21,30]
                             "chr1,21,+,40S19M,60,0;";      // DUP:TANDEM: [21,30]

        analyze_sa_tag(read_name, flag, chromosome, pos, mapq, test_cigar, seq, sa_tag, args, junctions_res);

        std::vector<Junction> junctions_expected_res =
        {                                                                     // +1 as we are 0 based but SAM is 1 based
            Junction{Breakend{"chr1", 20, strand::forward}, Breakend{"chr1", 29, strand::forward}, // DUP:TANDEM
                     "AAACGGCCCC"_dna5, 4, read_name}                                              // [21,30] duplicated
        };

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(),
                      junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i])
                        << "Junction " << i << " unequal\n"
                        << "Mate 1 equal: " << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                        << "\nMate 2 equal: " << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                        << "\n";
        }

        // For debugging use:
        // print_compare_junction_vectors(junctions_expected_res, junctions_res);
    }

    { // Example 3 - multiple tandem duplication on read and ref
        seqan3::debug_stream << "----------------------------------Third Example:----------------------------------\n";
        std::vector<Junction> junctions_res{};

        // pos  11  21           31                                                  41
        //      |   |            |                                                   |
        // chr1 -1->-----2------>-----2------>                                       -3->
        //      ||||||||||||||||||||||||||||||\\\\\\\\\\\\ \\\\\\\\\\\\ \\\\\\\\\\\\ ||||
        // read --->-DUP:TANDEM->-DUP:TANDEM->-DUP:TANDEM->-DUP:TANDEM->-DUP:TANDEM->--->
        //          |            |            |            |            |            |
        // read_pos 11           21           31           41           51           61

        // Primary alignment: chr1,11,+,30M39S,60,0;
        std::string read_name = "read";
        seqan3::sam_flag flag{0};
        std::string chromosome = "chr1";
        int32_t pos = 10; // Decrement by 1 because scanned position is already 0-based
        uint8_t mapq = 60;
        std::vector<seqan3::cigar> test_cigar = {{30, 'M'_cigar_operation}, {39, 'S'_cigar_operation}};
        seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5}; // 1 - 50

        // Supplementary alignments
        std::string sa_tag = /*"chr1,11,+,30M39S,60,0;"*/   // AS=30 -> primary;
                             "chr1,21,+,30S20M19S,60,0;"    // DUP:TANDEM: [21,40]
                             "chr1,31,+,50S19M,60,0;";      // DUP:TANDEM: [21,40]

        analyze_sa_tag(read_name, flag, chromosome, pos, mapq, test_cigar, seq, sa_tag, args, junctions_res);

        std::vector<Junction> junctions_expected_res =
        {                                                                     // +1 as we are 0 based but SAM is 1 based
            Junction{Breakend{"chr1", 20, strand::forward}, Breakend{"chr1", 39, strand::forward}, // DUP:TANDEM
                     ""_dna5, 5, read_name}                                                        // [21,40] duplicated
        };

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(),
                      junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i])
                        << "Junction " << i << " unequal\n"
                        << "Mate 1 equal: " << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                        << "\nMate 2 equal: " << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                        << "\n";
        }

        // For debugging use:
        // print_compare_junction_vectors(junctions_expected_res, junctions_res);
    }

    { // Example 4 - multiple tandem duplication on read and ref
        seqan3::debug_stream << "----------------------------------Fourth Example:----------------------------------\n";
        std::vector<Junction> junctions_res{};

        // pos  11  21           31           41                                     51
        //      |   |            |            |                                      |
        // chr1 -1->-----2------>-----2------>-----2------>                          -3->
        //      ||||||||||||||||||||||||||||||\\\\\\\\\\\\ \\\\\\\\\\\\ \\\\\\\\\\\\ ||||
        // read --->-DUP:TANDEM->-DUP:TANDEM->-DUP:TANDEM->-DUP:TANDEM->-DUP:TANDEM->--->
        //          |            |            |            |            |            |
        // read_pos 11           21           31           41           51           61

        // Primary alignment: chr1,11,+,40M29S,60,0;
        std::string read_name = "read";
        seqan3::sam_flag flag{0};
        std::string chromosome = "chr1";
        int32_t pos = 10; // Decrement by 1 because scanned position is already 0-based
        uint8_t mapq = 60;
        std::vector<seqan3::cigar> test_cigar = {{40, 'M'_cigar_operation}, {29, 'S'_cigar_operation}};
        seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"_dna5}; // 1 - 50

        // Supplementary alignments
        std::string sa_tag = /*"chr1,11,+,40M29S,60,0;"*/   // AS=30 -> primary;
                             "chr1,31,+,40S29M,60,0;";      // DUP:TANDEM: [21,50]

        analyze_sa_tag(read_name, flag, chromosome, pos, mapq, test_cigar, seq, sa_tag, args, junctions_res);

        // TODO (irallia 01.10.2021): How likely is this special case, do we need to solve this for CNV?
        // The algorithm currently finds only 2 duplications, since it does not perceive the first one in the reference
        // and thus interprets the other 2 as one longer one and returns a start position of 30 instead of 20.
        std::vector<Junction> junctions_expected_res =
        {                                                                     // +1 as we are 0 based but SAM is 1 based
            Junction{Breakend{"chr1", 30, strand::forward}, Breakend{"chr1", 49, strand::forward}, // DUP:TANDEM
                     "AAACGGCCCC"_dna5, 2, read_name}                                              // [21,50] duplicated
        };

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(),
                      junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i])
                        << "Junction " << i << " unequal\n"
                        << "Mate 1 equal: " << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                        << "\nMate 2 equal: " << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                        << "\n";
        }

        // For debugging use:
        // print_compare_junction_vectors(junctions_expected_res, junctions_res);
    }

    { // Example 5
        seqan3::debug_stream << "----------------------------------Fifth Example:----------------------------------\n";
        std::vector<Junction> junctions_res{};

        // chr1:                                        ??? GGGCTC TTCGGATCG GGCAGCATCAACGCT AAAC ???????????????? GGCCCC TGACAGGATA
        //                                                                  ||||||??||||||      (test_cigar: 60S14M26S)
        // read: GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCCGCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGACTGACAGGATA

        //                  | pos 107                                              | pos 157
        //           | pos 101         | pos 117      | pos 131             | pos 151
        // chr1: ??? GGGCTC TTCGGATCGG GCAGCATCAACGCT AAAC ???????????????? GGCCCC TGACAGGATA ???
        //           |||||| |||||||||| |||inversion|| |    ?|?|?|?|?|?|?|?|        ||||||||||   (SA-Tags: 6M94S & 16S10M74S & 60S14M26S & 90S10M)
        // read:     GGGCTC.TTCGGATCGG.TCGCAACTACGACG /    CGCATCAGCAGGCGAC        TGACAGGATA
        //                / \        / \             v AAAC GGCCCC moved to left                (SA-Tags: 40S4M56S & 44S6M50S)
        //            ATCGATCGAT     GGGGCCCCCATTTT  AAAC GGCCCC GCGATACGCG                     (insertion of GCGATACGCG?)
        //            ||||||||||                                                                (SA-Tag: 6S10M84S)
        // chr2:  ??? ATCGATCGAT ???                                                            (translocation of ATCGATCGAT chr2 -> chr1)
        //            | pos 101

        // Primary alignment: chr1,116,-,10S14M26S,60,0;
        std::string read_name = "read021";
        seqan3::sam_flag flag{16u}; // read reverse strand (0x10)
        std::string chromosome = "chr1";
        int32_t pos = 116;
        uint8_t mapq = 60;
        std::vector<seqan3::cigar> test_cigar = {{60, 'S'_cigar_operation},
                                                 {14, 'M'_cigar_operation},
                                                 {26, 'S'_cigar_operation}};
        seqan3::dna5_vector seq = {"GGGCTCATCGATCGATTTCGGATCGGGGGGCCCCCATTTTAAACGGCCCC"
                                   "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGACTGACAGGATA"_dna5};
        // Supplementary alignments
        std::string sa_tag = "chr1,101,+,6M94S,60,0;"
                             "chr2,101,+,6S10M84S,60,0;"
                             "chr1,107,+,16S10M74S,60,0;"
                             "chr1,117,-,60S14M26S,60,0;"
                             "chr1,131,+,40S4M56S,60,0;"
                             "chr1,151,+,44S6M50S,60,0;"
                             "chr1,157,+,90S10M,60,0;";

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
        std::vector<Junction> junctions_expected_res{Junction{new_breakend_1, new_breakend_2,
                                                              ""_dna5,
                                                              0, read_name},         // translocation
                                                     Junction{new_breakend_3, new_breakend_4,
                                                              ""_dna5,
                                                              0, read_name},         // translocation
                                                     Junction{new_breakend_5, new_breakend_6,
                                                              ""_dna5,
                                                              0, read_name},         // inversion
                                                     Junction{new_breakend_7, new_breakend_8,
                                                              ""_dna5,
                                                              0, read_name},         // inversion
                                                     Junction{new_breakend_9, new_breakend_10,
                                                              ""_dna5,
                                                              0, read_name},         // deletion
                                                     Junction{new_breakend_11, new_breakend_12,
                                                              "GCGATACGCGTCGCAACTACGACGCGCATCAGCAGGCGAC"_dna5,
                                                              0, read_name}          // insertion
                                                    };

        ASSERT_EQ(junctions_expected_res.size(), junctions_res.size());

        for (size_t i = 0; i < junctions_expected_res.size(); ++i)
        {
            EXPECT_EQ(junctions_expected_res[i].get_read_name(),
                      junctions_res[i].get_read_name()) << "Read names of junction " << i << " unequal";
            EXPECT_TRUE(junctions_expected_res[i] == junctions_res[i])
                        << "Junction " << i << " unequal\n"
                        << "Mate 1 equal: " << (junctions_expected_res[i].get_mate1() == junctions_res[i].get_mate1())
                        << "\nMate 2 equal: " << (junctions_expected_res[i].get_mate2() == junctions_res[i].get_mate2())
                        << "\n";
        }

        // For debugging use:
        // print_compare_junction_vectors(junctions_expected_res, junctions_res);
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
