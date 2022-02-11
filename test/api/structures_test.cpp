#include "api_test.hpp"

#include "structures/aligned_segment.hpp"
#include "structures/breakend.hpp"
#include "structures/cluster.hpp"
#include "structures/junction.hpp"
#include "variant_detection/method_enums.hpp"

/* tests for aligned_segment */

using seqan3::operator""_cigar_operation;
using seqan3::operator""_dna5;

TEST(structures, aligned_segment)
{
    // get_left_soft_clip()
    {
        // no soft clip with operation '='
        std::vector<seqan3::cigar> cig = {{1, '='_cigar_operation}};
        AlignedSegment aligned_segment = {strand::forward, "chr1", 100, 60, cig};
        EXPECT_EQ(aligned_segment.get_left_soft_clip(), 0);

        // no soft clip with operation 'I'
        cig = {{1, 'I'_cigar_operation}};
        aligned_segment = {strand::forward, "chr1", 100, 60, cig};
        EXPECT_EQ(aligned_segment.get_left_soft_clip(), 0);
    }

    // get_query_length()
    {
        // zero query length
        std::vector<seqan3::cigar> cig = {{1, 'D'_cigar_operation}};
        AlignedSegment aligned_segment = {strand::forward, "chr1", 100, 60, cig};
        EXPECT_EQ(aligned_segment.get_query_length(), 0);
    }
}

/* tests for method_enums */

TEST(structures, method_enums_detection_methods)
{
    std::unordered_map<std::string,
                       detection_methods> cigar_string_mapping = enumeration_names(detection_methods::cigar_string);
    EXPECT_EQ(detection_methods::cigar_string, cigar_string_mapping["cigar_string"]);
    EXPECT_EQ(detection_methods::cigar_string, cigar_string_mapping["0"]);

    std::unordered_map<std::string,
                       detection_methods> split_read_mapping = enumeration_names(detection_methods::split_read);
    EXPECT_EQ(detection_methods::split_read, split_read_mapping["split_read"]);
    EXPECT_EQ(detection_methods::split_read, split_read_mapping["1"]);

    std::unordered_map<std::string,
                       detection_methods> read_pairs_mapping = enumeration_names(detection_methods::read_pairs);
    EXPECT_EQ(detection_methods::read_pairs, read_pairs_mapping["read_pairs"]);
    EXPECT_EQ(detection_methods::read_pairs, read_pairs_mapping["2"]);

    std::unordered_map<std::string,
                       detection_methods> read_depth_mapping = enumeration_names(detection_methods::read_depth);
    EXPECT_EQ(detection_methods::read_depth, read_depth_mapping["read_depth"]);
    EXPECT_EQ(detection_methods::read_depth, read_depth_mapping["3"]);
}

TEST(structures, method_enums_clustering_methods)
{
    std::unordered_map<std::string, clustering_methods> simple_clustering_mapping
                                        = enumeration_names(clustering_methods::simple_clustering);
    EXPECT_EQ(clustering_methods::simple_clustering, simple_clustering_mapping["simple_clustering"]);
    EXPECT_EQ(clustering_methods::simple_clustering, simple_clustering_mapping["0"]);

    std::unordered_map<std::string, clustering_methods> hierarchical_clustering_mapping
                                        = enumeration_names(clustering_methods::hierarchical_clustering);
    EXPECT_EQ(clustering_methods::hierarchical_clustering, hierarchical_clustering_mapping["hierarchical_clustering"]);
    EXPECT_EQ(clustering_methods::hierarchical_clustering, hierarchical_clustering_mapping["1"]);

    std::unordered_map<std::string, clustering_methods> self_balancing_binary_tree_mapping
                                        = enumeration_names(clustering_methods::self_balancing_binary_tree);
    EXPECT_EQ(clustering_methods::self_balancing_binary_tree,
              self_balancing_binary_tree_mapping["self_balancing_binary_tree"]);
    EXPECT_EQ(clustering_methods::self_balancing_binary_tree, self_balancing_binary_tree_mapping["2"]);

    std::unordered_map<std::string, clustering_methods> candidate_selection_based_on_voting_mapping
                                        = enumeration_names(clustering_methods::candidate_selection_based_on_voting);
    EXPECT_EQ(clustering_methods::candidate_selection_based_on_voting,
              candidate_selection_based_on_voting_mapping["candidate_selection_based_on_voting"]);
    EXPECT_EQ(clustering_methods::candidate_selection_based_on_voting,
              candidate_selection_based_on_voting_mapping["3"]);
}

TEST(structures, method_enums_refinement_methods)
{
    std::unordered_map<std::string, refinement_methods> no_refinement_mapping
                                        = enumeration_names(refinement_methods::no_refinement);
    EXPECT_EQ(refinement_methods::no_refinement, no_refinement_mapping["no_refinement"]);
    EXPECT_EQ(refinement_methods::no_refinement, no_refinement_mapping["0"]);

    std::unordered_map<std::string, refinement_methods> sViper_refinement_method_mapping
                                        = enumeration_names(refinement_methods::sViper_refinement_method);
    EXPECT_EQ(refinement_methods::sViper_refinement_method,
              sViper_refinement_method_mapping["sViper_refinement_method"]);
    EXPECT_EQ(refinement_methods::sViper_refinement_method, sViper_refinement_method_mapping["1"]);

    std::unordered_map<std::string, refinement_methods> sVirl_refinement_method_mapping
                                        = enumeration_names(refinement_methods::sVirl_refinement_method);
    EXPECT_EQ(refinement_methods::sVirl_refinement_method, sVirl_refinement_method_mapping["sVirl_refinement_method"]);
    EXPECT_EQ(refinement_methods::sVirl_refinement_method, sVirl_refinement_method_mapping["2"]);
}

/* tests for junctions */

TEST(structures, junctions_breakend_flip_orientation)
{
    Breakend forward_breakend{"chr1", 42, strand::forward};
    Breakend reverse_breakend{"chr1", 42, strand::reverse};

    EXPECT_NE(forward_breakend, reverse_breakend);
    forward_breakend.flip_orientation();
    EXPECT_EQ(forward_breakend, reverse_breakend); // both are reverse now
    reverse_breakend.flip_orientation();
    EXPECT_NE(forward_breakend, reverse_breakend);
    forward_breakend.flip_orientation();
    EXPECT_EQ(forward_breakend, reverse_breakend); // both are forward now
}

/* tests for clusters */

TEST(structures, clusters_get_common_tandem_dup_count)
{
    Junction junction_1{Breakend{"chr1", 123456, strand::forward},
                        Breakend{"chr1", 234567, strand::forward},
                        ""_dna5, 0, "read_1"};
    Junction junction_2{Breakend{"chr1", 123456, strand::forward},
                        Breakend{"chr1", 234567, strand::forward},
                        ""_dna5, 1, "read_1"};
    Junction junction_3{Breakend{"chr1", 123456, strand::forward},
                        Breakend{"chr1", 234567, strand::forward},
                        ""_dna5, 2, "read_1"};
    Junction junction_4{Breakend{"chr1", 123456, strand::forward},
                        Breakend{"chr1", 234567, strand::forward},
                        ""_dna5, 7, "read_1"};

    // Tandem duplication with 2 copies:
    Cluster cluster_1{{junction_3, junction_3, junction_3}};
    EXPECT_EQ(cluster_1.get_common_tandem_dup_count(), 2);

    // Tandem duplication with different amount of copies:
    Cluster cluster_2{{junction_2, junction_3, junction_4}};
    EXPECT_EQ(cluster_2.get_common_tandem_dup_count(), 3); // (1+2+7)/3=3,333

    // More than two thirds of the junctions have a 0 tandem_dup_count, than its probably no tandem duplication.
    Cluster cluster_3{{junction_1, junction_1, junction_1, junction_2}};
    EXPECT_EQ(cluster_3.get_common_tandem_dup_count(), 0);

    // Two thirds of the junctions have a 0 tandem_dup_count, than its probably a tandem duplication.
    Cluster cluster_4{{junction_1, junction_1, junction_2}};
    EXPECT_EQ(cluster_4.get_common_tandem_dup_count(), 0);

    // Less than two thirds of the junctions have a 0 tandem_dup_count, than its probably a tandem duplication.
    Cluster cluster_5{{junction_1, junction_2}};
    EXPECT_EQ(cluster_5.get_common_tandem_dup_count(), 1);
}
