#include <gtest/gtest.h>

#include "structures/breakend.hpp"
#include "variant_detection/method_enums.hpp"

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

TEST(structures, breakend_flip_orientation)
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
