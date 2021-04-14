#include <gtest/gtest.h>

#include "modules/clustering/simple_clustering_method.hpp"  // for the simple clustering method
#include "structures/cluster.hpp"                           // for class Cluster

/* -------- clustering methods tests -------- */

TEST(junction_detection, clustering_method_simple)
{
    testing::internal::CaptureStdout();

    std::vector<Junction> junctions{};
    std::vector<Cluster> resulting_clusters{};
    simple_clustering_method(junctions, resulting_clusters);

    std::vector<Cluster> expected_clusters{};
    // TODO (irallia): We need to implement operator== for Clusters.
    // EXPECT_EQ(expected_clusters, resulting_clusters);
}
