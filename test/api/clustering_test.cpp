#include <gtest/gtest.h>

#include "modules/clustering/simple_clustering_method.hpp"          // for the simple clustering method
#include "modules/clustering/hierarchical_clustering_method.hpp"    // for the hierarchical clustering method
#include "structures/cluster.hpp"                                   // for class Cluster

using seqan3::operator""_dna5;

/* -------- clustering methods tests -------- */

std::string const chrom1 = "chr1";
int32_t const chrom1_position1 = 12323443;
int32_t const chrom1_position2 = 94734377;
int32_t const chrom1_position3 = 112323345;
std::string const chrom2 = "chr2";
int32_t const chrom2_position1 = 234432;
std::string const read_name_1 = "m2257/8161/CCS";
std::string const read_name_2 = "m41327/11677/CCS";
std::string const read_name_3 = "m21263/13017/CCS";
std::string const read_name_4 = "m38637/7161/CCS";
std::string const read_name_5 = "m23412/9534/CCS";
std::string const read_name_6 = "m1245/5634/CCS";
std::string const read_name_7 = "m8765/9765/CCS";
std::string const read_name_8 = "m13456/11102/CCS";

std::vector<Junction> prepare_input_junctions()
{
    std::vector<Junction> input_junctions
    {
        Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                 Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, read_name_1}, //cluster 1
        Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                 Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, read_name_2}, //cluster 1
        Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                 Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, read_name_3}, //cluster 1
        Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                 Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, read_name_4}, //cluster 2
        Junction{Breakend{chrom1, chrom1_position1 + 92, strand::forward},
                 Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, read_name_5}, //cluster 3
        Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                 Breakend{chrom2, chrom1_position3 + 8, strand::reverse}, ""_dna5, read_name_6}, //cluster 4
        Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                 Breakend{chrom2, chrom1_position3 - 1, strand::reverse}, ""_dna5, read_name_7}, //cluster 4
        Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                 Breakend{chrom2, chrom1_position3 + 2, strand::reverse}, ""_dna5, read_name_8} //cluster 4
    };

    std::sort(input_junctions.begin(), input_junctions.end());
    return input_junctions;
}

TEST(clustering, clustering_method_simple)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> resulting_clusters{};
    resulting_clusters = simple_clustering_method(input_junctions);

    // Each junction in separate cluster
    std::vector<Cluster> expected_clusters
    {
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, read_name_1}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, read_name_2}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, read_name_3},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, read_name_4},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 92, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, read_name_5},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 8, strand::reverse}, ""_dna5, read_name_6}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                             Breakend{chrom2, chrom1_position3 - 1, strand::reverse}, ""_dna5, read_name_7}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 2, strand::reverse}, ""_dna5, read_name_8}
        }}
    };
    std::sort(expected_clusters.begin(), expected_clusters.end());

    ASSERT_EQ(expected_clusters.size(), resulting_clusters.size());

    for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
    {
        EXPECT_TRUE(expected_clusters[cluster_index] == resulting_clusters[cluster_index]) << "Cluster " << cluster_index << " unequal";
    }
}

TEST(clustering, partitioning)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<std::vector<Junction>> partitions = partition_junctions(input_junctions);

    std::vector<std::vector<Junction>> expected_partitions
    {
        {
            Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, read_name_1}, //cluster 1
            Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                     Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, read_name_2}, //cluster 1
            Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, read_name_3}, //cluster 1
        },
        {
            Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                     Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, read_name_4}, //cluster 2
        },
        {
            Junction{Breakend{chrom1, chrom1_position1 + 92, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, read_name_5}, //cluster 3
        },
        {
            Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                     Breakend{chrom2, chrom1_position3 + 8, strand::reverse}, ""_dna5, read_name_6}, //cluster 4
            Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                     Breakend{chrom2, chrom1_position3 - 1, strand::reverse}, ""_dna5, read_name_7}, //cluster 4
            Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                     Breakend{chrom2, chrom1_position3 + 2, strand::reverse}, ""_dna5, read_name_8} //cluster 4
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


TEST(clustering, strict_clustering)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, 0);

    // Each junction in separate cluster
    std::vector<Cluster> expected_clusters
    {
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, read_name_1}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, read_name_2}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, read_name_3},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, read_name_4},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 92, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, read_name_5},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 8, strand::reverse}, ""_dna5, read_name_6}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                             Breakend{chrom2, chrom1_position3 - 1, strand::reverse}, ""_dna5, read_name_7}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 2, strand::reverse}, ""_dna5, read_name_8}
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


TEST(clustering, clustering_10)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, 10);

    // Distance matrix for junctions from reads 1-3 and 6-8
    //      1   2   3
    //  1       18  21
    //  2           11
    //
    //      6   7   8
    //  6       14  14
    //  7           6

    // Only junctions from reads 7 and 8 have a distance < 10 and cluster together
    std::vector<Cluster> expected_clusters
    {
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, read_name_1}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, read_name_2}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, read_name_3},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, read_name_4},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 92, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, read_name_5},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 8, strand::reverse}, ""_dna5, read_name_6}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                             Breakend{chrom2, chrom1_position3 - 1, strand::reverse}, ""_dna5, read_name_7},
                    Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 2, strand::reverse}, ""_dna5, read_name_8}
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


TEST(clustering, clustering_15)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, 15);

    // Distance matrix for junctions from reads 1-3 and 6-8
    //      1   2   3
    //  1       18  21
    //  2           11
    //
    //      6   7   8
    //  6       14  14
    //  7           6

    // Junctions from reads 6-8 and 2-3 have a distance < 15 and cluster together
    std::vector<Cluster> expected_clusters
    {
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, read_name_1}
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, read_name_2},
                    Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, read_name_3},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, read_name_4},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 92, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, read_name_5},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 8, strand::reverse}, ""_dna5, read_name_6},
                    Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                             Breakend{chrom2, chrom1_position3 - 1, strand::reverse}, ""_dna5, read_name_7},
                    Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 2, strand::reverse}, ""_dna5, read_name_8}
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


TEST(clustering, clustering_25)
{
    std::vector<Junction> input_junctions = prepare_input_junctions();
    std::vector<Cluster> clusters = hierarchical_clustering_method(input_junctions, 25);

    // Distance matrix for junctions from reads 1-3 and 6-8
    //      1   2   3
    //  1       18  21
    //  2           11
    //
    //      6   7   8
    //  6       14  14
    //  7           6

    // Junctions from reads 6-8 and 1-3 have a distance < 25 and cluster together
    std::vector<Cluster> expected_clusters
    {
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, read_name_1},
                    Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, read_name_2},
                    Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, read_name_3},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                             Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, read_name_4},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position1 + 92, strand::forward},
                             Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, read_name_5},
        }},
        Cluster{{   Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 8, strand::reverse}, ""_dna5, read_name_6},
                    Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                             Breakend{chrom2, chrom1_position3 - 1, strand::reverse}, ""_dna5, read_name_7},
                    Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                             Breakend{chrom2, chrom1_position3 + 2, strand::reverse}, ""_dna5, read_name_8}
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
