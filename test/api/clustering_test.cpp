#include <gtest/gtest.h>

#include "modules/clustering/simple_clustering_method.hpp"          // for the simple clustering method
#include "modules/clustering/hierarchical_clustering_method.hpp"    // for the hierarchical clustering method
#include "structures/cluster.hpp"                                   // for class Cluster

using seqan3::operator""_dna5;

/* -------- clustering methods tests -------- */

TEST(clustering, clustering_method_simple)
{
    testing::internal::CaptureStdout();

    std::vector<Junction> junctions{};
    std::vector<Cluster> resulting_clusters{};
    simple_clustering_method(junctions, resulting_clusters);

    std::vector<Cluster> expected_clusters{};
    // TODO (irallia): We need to implement operator== for Clusters.
    // EXPECT_EQ(expected_clusters, resulting_clusters);
}

TEST(clustering, partitioning)
{
    const std::string chrom1 = "chr1";
    const int32_t chrom1_position1 = 12323443;
    const int32_t chrom1_position2 = 94734377;
    const int32_t chrom1_position3 = 112323345;
    const std::string chrom2 = "chr2";
    const int32_t chrom2_position1 = 234432;
    const std::string & read_name_1 = "m2257/8161/CCS";
    const std::string & read_name_2 = "m41327/11677/CCS";
    const std::string & read_name_3 = "m21263/13017/CCS";
    const std::string & read_name_4 = "m38637/7161/CCS";
    const std::string & read_name_5 = "m23412/9534/CCS";
    const std::string & read_name_6 = "m1245/5634/CCS";
    const std::string & read_name_7 = "m8765/9765/CCS";
    const std::string & read_name_8 = "m13456/11102/CCS";

    std::vector<Junction> input_junctions
    {   Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
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

    std::vector<std::vector<Junction>> clusters = partition_junctions(input_junctions);

    std::vector<std::vector<Junction>> expected_clusters
    {   
        {   Junction{Breakend{chrom1, chrom1_position1 - 5, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 8, strand::forward}, ""_dna5, read_name_1}, //cluster 1
            Junction{Breakend{chrom1, chrom1_position1 + 2, strand::forward},
                     Breakend{chrom2, chrom2_position1 - 3, strand::forward}, ""_dna5, read_name_2}, //cluster 1
            Junction{Breakend{chrom1, chrom1_position1 + 9, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 1, strand::forward}, ""_dna5, read_name_3}, //cluster 1
        },
        {   Junction{Breakend{chrom1, chrom1_position1 + 5, strand::forward},
                     Breakend{chrom2, chrom2_position1 - 1, strand::reverse}, ""_dna5, read_name_4}, //cluster 2
        },
        {   Junction{Breakend{chrom1, chrom1_position1 + 92, strand::forward},
                     Breakend{chrom2, chrom2_position1 + 3, strand::forward}, ""_dna5, read_name_5}, //cluster 3
        },
        {   Junction{Breakend{chrom1, chrom1_position2 - 2, strand::forward},
                     Breakend{chrom2, chrom1_position3 + 8, strand::reverse}, ""_dna5, read_name_6}, //cluster 4
            Junction{Breakend{chrom1, chrom1_position2 + 3, strand::forward},
                     Breakend{chrom2, chrom1_position3 - 1, strand::reverse}, ""_dna5, read_name_7}, //cluster 4
            Junction{Breakend{chrom1, chrom1_position2 + 6, strand::forward},
                     Breakend{chrom2, chrom1_position3 + 2, strand::reverse}, ""_dna5, read_name_8} //cluster 4
        }
    };

    EXPECT_EQ(expected_clusters.size(), clusters.size());

    if (expected_clusters.size() == clusters.size())
    {
        for (size_t cluster_index = 0; cluster_index < expected_clusters.size(); ++cluster_index)
        {
            EXPECT_EQ(expected_clusters[cluster_index].size(), clusters[cluster_index].size());
            if (expected_clusters[cluster_index].size() == clusters[cluster_index].size())
            {
                for (size_t junction_index = 0; junction_index < expected_clusters[cluster_index].size(); ++junction_index)
                {
                    EXPECT_TRUE(expected_clusters[cluster_index][junction_index] == clusters[cluster_index][junction_index]);
                }
            }
        }
    }
}
