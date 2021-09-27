#include <gtest/gtest.h>

#include "iGenVar.hpp"                                              // for global variable gVerbose
#include "modules/clustering/simple_clustering_method.hpp"          // for the simple clustering method
#include "modules/clustering/hierarchical_clustering_method.hpp"    // for the hierarchical clustering method
#include "structures/cluster.hpp"                                   // for class Cluster

using seqan3::operator""_dna5;

/* -------- clustering methods tests -------- */

std::string const chrom1 = "chr1";
int32_t const chrom1_position1 = 12323443;
int32_t const chrom1_position2 = 94734377;
std::string const chrom2 = "chr2";
int32_t const chrom2_position1 = 234432;
int32_t const chrom2_position2 = 112323345;
size_t tandem_dup_count = 0;
std::string const read_name_1 = "m2257/8161/CCS";
std::string const read_name_2 = "m41327/11677/CCS";
std::string const read_name_3 = "m21263/13017/CCS";
std::string const read_name_4 = "m38637/7161/CCS";
std::string const read_name_5 = "m23412/9534/CCS";
std::string const read_name_6 = "m1245/5634/CCS";
std::string const read_name_7 = "m8765/9765/CCS";
std::string const read_name_8 = "m13456/11102/CCS";

constexpr int32_t default_partition_max_distance = 1000;

std::vector<Junction> prepare_input_junctions()
{
    std::vector<Junction> input_junctions
    {
        Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                 Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, tandem_dup_count, read_name_1},
        Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                 Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, tandem_dup_count, read_name_2},
        Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                 Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, tandem_dup_count, read_name_3},
        Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                 Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, tandem_dup_count, read_name_4},
        Junction{Breakend{chrom1, chrom1_position1 + 1245, strand::forward},
                 Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, tandem_dup_count, read_name_5},
        Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                 Breakend{chrom2, chrom2_position2 + 8, strand::reverse}, ""_dna5, tandem_dup_count, read_name_6},
        Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                 Breakend{chrom2, chrom2_position2 - 1, strand::reverse}, ""_dna5, tandem_dup_count, read_name_7},
        Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                 Breakend{chrom2, chrom2_position2 + 2, strand::reverse}, ""_dna5, tandem_dup_count, read_name_8}
    };

    std::sort(input_junctions.begin(), input_junctions.end());
    return input_junctions;
}

TEST(simple_clustering, all_separate)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> resulting_clusters{};
    resulting_clusters = simple_clustering_method(input_junctions);

    // Each junction in separate cluster
    std::vector<Cluster> expected_clusters
    {
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 8, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_1}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_2}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 1, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_3},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_4},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 1245, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_5},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 8, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_6}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                          Breakend{chrom2, chrom2_position2 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_7}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 2, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_8}
        }}
    };
    std::sort(expected_clusters.begin(), expected_clusters.end());

    ASSERT_EQ(expected_clusters.size(), resulting_clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
    {
        EXPECT_TRUE(expected_clusters[cluster_index] == resulting_clusters[cluster_index]) << "Cluster "
                                                                                           << cluster_index
                                                                                           << " unequal";
    }
}

TEST(simple_clustering, clustered)
{
    std::vector<Junction> input_junctions
    {
        Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                 Breakend{chrom2, chrom2_position1, strand::forward}, ""_dna5, tandem_dup_count, read_name_1},
        Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                 Breakend{chrom2, chrom2_position1, strand::forward}, ""_dna5, tandem_dup_count, read_name_2}
    };
    std::vector<Cluster> resulting_clusters{};
    resulting_clusters = simple_clustering_method(input_junctions);

    // Both junctions in the same cluster
    std::vector<Cluster> expected_clusters
    {
        Cluster{{Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                          Breakend{chrom2, chrom2_position1, strand::forward}, ""_dna5, tandem_dup_count, read_name_1},
                 Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                          Breakend{chrom2, chrom2_position1, strand::forward}, ""_dna5, tandem_dup_count, read_name_2}
        }}
    };
    std::sort(expected_clusters.begin(), expected_clusters.end());

    ASSERT_EQ(1u, resulting_clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
    {
        EXPECT_TRUE(expected_clusters[cluster_index] == resulting_clusters[cluster_index]) << "Cluster "
                                                                                           << cluster_index
                                                                                           << " unequal";
    }
}

TEST(simple_clustering, empty_junction_vector)
{
    testing::internal::CaptureStdout();
    testing::internal::CaptureStderr();

    std::vector<Junction> input_junctions{};
    std::vector<Cluster> resulting_clusters{};
    resulting_clusters = simple_clustering_method(input_junctions);

    ASSERT_EQ(0u, resulting_clusters.size());

    std::string result_out = testing::internal::GetCapturedStdout();
    EXPECT_EQ("", result_out);
    std::string result_err = testing::internal::GetCapturedStderr();
    EXPECT_EQ("No junctions found...\n", result_err);
}

TEST(hierarchical_clustering, partitioning)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<std::vector<Junction>> partitions = partition_junctions(input_junctions, default_partition_max_distance);

    std::vector<std::vector<Junction>> expected_partitions
    {
        {
            Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, tandem_dup_count, read_name_1}, //cluster 1
            Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                     Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, tandem_dup_count, read_name_2}, //cluster 1
            Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, tandem_dup_count, read_name_3}, //cluster 1
        },
        {
            Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                     Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, tandem_dup_count, read_name_4}, //cluster 2
        },
        {
            Junction{Breakend{chrom1, chrom1_position1 + 1245, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, tandem_dup_count, read_name_5}, //cluster 3
        },
        {
            Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                     Breakend{chrom2, chrom2_position2 + 8, strand::reverse}, ""_dna5, tandem_dup_count, read_name_6}, //cluster 4
            Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                     Breakend{chrom2, chrom2_position2 - 1, strand::reverse}, ""_dna5, tandem_dup_count, read_name_7}, //cluster 4
            Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                     Breakend{chrom2, chrom2_position2 + 2, strand::reverse}, ""_dna5, tandem_dup_count, read_name_8}  //cluster 4
        }
    };

    ASSERT_EQ(expected_partitions.size(), partitions.size());

    for (size_t partition_index = 0; partition_index < expected_partitions.size(); ++partition_index)
    {
        ASSERT_EQ(expected_partitions[partition_index].size(), partitions[partition_index].size());
        for (size_t junction_index = 0; junction_index < expected_partitions[partition_index].size(); ++junction_index)
        {
            EXPECT_TRUE(expected_partitions[partition_index][junction_index] == partitions[partition_index][junction_index]);
        }
    }
}

TEST(hierarchical_clustering, strict_clustering)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, default_partition_max_distance, 0);

    // Each junction in separate cluster
    std::vector<Cluster> expected_clusters
    {
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 8, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_1}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_2}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 1, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_3},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_4},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 1245, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_5},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 8, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_6}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                          Breakend{chrom2, chrom2_position2 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_7}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 2, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_8}
        }}
    };
    std::sort(expected_clusters.begin(), expected_clusters.end());

    ASSERT_EQ(expected_clusters.size(), clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
    {
        ASSERT_EQ(expected_clusters[cluster_index].get_cluster_size(),
                  clusters[cluster_index].get_cluster_size())
                  << "Cluster " << cluster_index << " of unexpected size";

        for (size_t junction_index = 0;
             junction_index < expected_clusters[cluster_index].get_cluster_size();
             ++junction_index)
        {
            EXPECT_EQ(expected_clusters[cluster_index].get_members()[junction_index],
                      clusters[cluster_index].get_members()[junction_index])
                      << "Junction " << junction_index << " in cluster "
                      << cluster_index << " different than expected";
        }
    }
}

TEST(hierarchical_clustering, clustering_10)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, default_partition_max_distance, 0.01);

    // Position distance matrix for junctions from reads 1-3 and 6-8 (size distance is 0 for all pairs)
    //      1   2   3
    //  1       18  21
    //  2           11
    //
    //      6   7   8
    //  6       14  14
    //  7           6

    // Only junctions from reads 7 and 8 have a distance < 0.01 and cluster together
    std::vector<Cluster> expected_clusters
    {
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 8, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_1}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_2}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 1, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_3},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_4},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 1245, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_5},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 8, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_6}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                          Breakend{chrom2, chrom2_position2 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_7},
                 Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 2, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_8}
        }}
    };
    std::sort(expected_clusters.begin(), expected_clusters.end());

    ASSERT_EQ(expected_clusters.size(), clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
    {
        ASSERT_EQ(expected_clusters[cluster_index].get_cluster_size(),
                  clusters[cluster_index].get_cluster_size())
                  << "Cluster " << cluster_index << " of unexpected size";

        for (size_t junction_index = 0;
             junction_index < expected_clusters[cluster_index].get_cluster_size();
             ++junction_index)
        {
            EXPECT_EQ(expected_clusters[cluster_index].get_members()[junction_index],
                      clusters[cluster_index].get_members()[junction_index])
                      << "Junction " << junction_index << " in cluster "
                      << cluster_index << " different than expected";
        }
    }
}

TEST(hierarchical_clustering, clustering_15)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, default_partition_max_distance, 0.015);

    // Position distance matrix for junctions from reads 1-3 and 6-8 (size distance is 0 for all pairs)
    //      1   2   3
    //  1       18  21
    //  2           11
    //
    //      6   7   8
    //  6       14  14
    //  7           6

    // Junctions from reads 6-8 and 2-3 have a distance < 0.015 and cluster together
    std::vector<Cluster> expected_clusters
    {
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 8, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_1}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_2},
                 Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 1, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_3},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_4},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 1245, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_5},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 8, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_6},
                 Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                          Breakend{chrom2, chrom2_position2 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_7},
                 Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 2, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_8}
        }}
    };
    std::sort(expected_clusters.begin(), expected_clusters.end());

    ASSERT_EQ(expected_clusters.size(), clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
    {
        ASSERT_EQ(expected_clusters[cluster_index].get_cluster_size(),
                  clusters[cluster_index].get_cluster_size())
                  << "Cluster " << cluster_index << " of unexpected size";

        for (size_t junction_index = 0;
             junction_index < expected_clusters[cluster_index].get_cluster_size();
             ++junction_index)
        {
            EXPECT_EQ(expected_clusters[cluster_index].get_members()[junction_index],
                      clusters[cluster_index].get_members()[junction_index])
                      << "Junction " << junction_index << " in cluster "
                      << cluster_index << " different than expected";
        }
    }
}

TEST(hierarchical_clustering, clustering_25)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, default_partition_max_distance, 0.025);

    // Position distance matrix for junctions from reads 1-3 and 6-8 (size distance is 0 for all pairs)
    //      1   2   3
    //  1       18  21
    //  2           11
    //
    //      6   7   8
    //  6       14  14
    //  7           6

    // Junctions from reads 6-8 and 1-3 have a distance < 0.025 and cluster together
    std::vector<Cluster> expected_clusters
    {
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 8, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_1},
                 Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_2},
                 Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 1, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_3},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                          Breakend{chrom2, chrom2_position1 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_4},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1 + 1245, strand::forward},
                          Breakend{chrom2, chrom2_position1 + 3, strand::forward},
                          ""_dna5, tandem_dup_count, read_name_5},
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 8, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_6},
                 Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                          Breakend{chrom2, chrom2_position2 - 1, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_7},
                 Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                          Breakend{chrom2, chrom2_position2 + 2, strand::reverse},
                          ""_dna5, tandem_dup_count, read_name_8}
        }}
    };
    std::sort(expected_clusters.begin(), expected_clusters.end());

    ASSERT_EQ(expected_clusters.size(), clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
    {
        ASSERT_EQ(expected_clusters[cluster_index].get_cluster_size(),
                  clusters[cluster_index].get_cluster_size())
                  << "Cluster " << cluster_index << " of unexpected size";

        for (size_t junction_index = 0;
             junction_index < expected_clusters[cluster_index].get_cluster_size();
             ++junction_index)
        {
            EXPECT_EQ(expected_clusters[cluster_index].get_members()[junction_index],
                      clusters[cluster_index].get_members()[junction_index])
                      << "Junction " << junction_index << " in cluster "
                      << cluster_index << " different than expected";
        }
    }
}

TEST(hierarchical_clustering, subsampling)
{
    gVerbose = true;

    std::vector<Junction> input_junctions;
    for (int32_t i = 0; i < 300; ++i)
    {
        input_junctions.emplace_back(Breakend{chrom1, chrom1_position1 + i, strand::forward},
                                     Breakend{chrom2, chrom2_position1 + i, strand::forward},
                                     ""_dna5,
                                     tandem_dup_count,
                                     read_name_1);
    }
    std::sort(input_junctions.begin(), input_junctions.end());

    testing::internal::CaptureStderr();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, default_partition_max_distance, 0);

    size_t num_junctions = 0;
    for (Cluster const & cluster : clusters)
    {
        num_junctions += cluster.get_cluster_size();
    }
    EXPECT_EQ(num_junctions, 200u);

    std::string const expected_err
    {
        "A partition exceeds the maximum size (300>200) and has to be subsampled. "
        "Representative partition member:\n"
        "[chr1\t12323443\tForward] -> [chr2\t234432\tForward]\n"
    };
    std::string result_err = testing::internal::GetCapturedStderr();
    EXPECT_EQ(expected_err, result_err);
}

TEST(hierarchical_clustering, cluster_tandem_dup_count)
{
    std::vector<Junction> input_junctions
    {
        Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                 Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 0, read_name_1},
        Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                 Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 0, read_name_2},
        Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                 Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 1, read_name_3},
        Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                 Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 2, read_name_4},
        Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                 Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 2, read_name_5}

    };
    std::vector<Cluster> resulting_clusters{};

    std::vector<Cluster> expected_clusters
    {
        Cluster{{Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                          Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 0, read_name_1},
                 Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                          Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 0, read_name_2}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                          Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 1, read_name_3}
        }},
        Cluster{{Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                          Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 2, read_name_4},
                 Junction{Breakend{chrom1, chrom1_position1, strand::forward},
                          Breakend{chrom1, chrom1_position2, strand::forward}, ""_dna5, 2, read_name_5}
        }}
    };
    std::sort(expected_clusters.begin(), expected_clusters.end());

    resulting_clusters = simple_clustering_method(input_junctions);
    ASSERT_EQ(expected_clusters.size(), resulting_clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
    {
        EXPECT_TRUE(expected_clusters[cluster_index] == resulting_clusters[cluster_index]) << "Cluster "
                                                                                           << cluster_index
                                                                                           << " unequal";
    }

    std::vector<Cluster> expected_clusters_2 { input_junctions };

    resulting_clusters = hierarchical_clustering_method(input_junctions, default_partition_max_distance, 0);   // clustering_cutoff = 0
    ASSERT_EQ(5u, resulting_clusters.size());

    resulting_clusters = hierarchical_clustering_method(input_junctions, default_partition_max_distance, 1);   // clustering_cutoff = 1
    ASSERT_EQ(expected_clusters_2.size(), resulting_clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters_2.size(); ++cluster_index)
    {
        EXPECT_TRUE(expected_clusters_2[cluster_index] == resulting_clusters[cluster_index]) << "Cluster "
                                                                                             << cluster_index
                                                                                             << " unequal";
    }

    resulting_clusters = hierarchical_clustering_method(input_junctions, default_partition_max_distance, 10);  // clustering_cutoff = 10 (default value)
    ASSERT_EQ(1u, resulting_clusters.size());
}
