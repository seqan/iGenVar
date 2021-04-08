#include "modules/clustering/hierarchical_clustering_method.hpp"

#include <limits>                                                 // for infinity

#include "fastcluster.h"                                          // for hclust_fast


std::vector<std::vector<Junction>> partition_junctions(std::vector<Junction> const & junctions)
{
    // Partition based on mate 1
    std::vector<Junction> current_partition{};
    std::vector<std::vector<Junction>> current_partition_splitted;
    std::vector<std::vector<Junction>> final_partitions{};

    for (Junction junction : junctions)
    {
        if (current_partition.empty())
        {
            current_partition.push_back(junction);
        }
        else
        {
            if (junction.get_mate1().seq_name != current_partition.back().get_mate1().seq_name ||
                junction.get_mate1().orientation != current_partition.back().get_mate1().orientation ||
                abs(junction.get_mate1().position - current_partition.back().get_mate1().position) > 50)
            {
                // Partition based on mate 2
                std::sort(current_partition.begin(), current_partition.end(), [](Junction a, Junction b) {
                    return a.get_mate2() < b.get_mate2();
                });
                current_partition_splitted = split_partition_based_on_mate2(current_partition);
                for (std::vector<Junction> partition : current_partition_splitted)
                {
                    final_partitions.push_back(partition);
                }
                current_partition.clear();
            }
            current_partition.push_back(junction);
        }        
    }
    if (!current_partition.empty())
    {
        std::sort(current_partition.begin(), current_partition.end(), [](Junction a, Junction b) {
            return a.get_mate2() < b.get_mate2();
        });
        current_partition_splitted = split_partition_based_on_mate2(current_partition);
        for (std::vector<Junction> partition : current_partition_splitted)
        {
            final_partitions.push_back(partition);
        }
    }
    return final_partitions;
}


// Partition based on mate2
std::vector<std::vector<Junction>> split_partition_based_on_mate2(std::vector<Junction> const & partition)
{
    std::vector<Junction> current_partition{};
    std::vector<std::vector<Junction>> splitted_partition{};

    for (Junction junction : partition)
    {
        if (current_partition.empty())
        {
            current_partition.push_back(junction);
        }
        else
        {
            if (junction.get_mate2().seq_name != current_partition.back().get_mate2().seq_name ||
                junction.get_mate2().orientation != current_partition.back().get_mate2().orientation ||
                abs(junction.get_mate2().position - current_partition.back().get_mate2().position) > 50)
            {
                std::sort(current_partition.begin(), current_partition.end());
                splitted_partition.push_back(current_partition);
                current_partition.clear();
            }
            current_partition.push_back(junction);
        }
    }
    if (!current_partition.empty())
    {
        std::sort(current_partition.begin(), current_partition.end());
        splitted_partition.push_back(current_partition);
    }    
    return splitted_partition;
}

int junction_distance(Junction const & lhs, Junction const & rhs)
{
    if ((lhs.get_mate1().seq_name == rhs.get_mate1().seq_name) &&
        (lhs.get_mate1().orientation == rhs.get_mate1().orientation) &&
        (lhs.get_mate2().seq_name == rhs.get_mate2().seq_name) &&
        (lhs.get_mate2().orientation == rhs.get_mate2().orientation))
    {
        return (std::abs(lhs.get_mate1().position - rhs.get_mate1().position) +
                std::abs(lhs.get_mate2().position - rhs.get_mate2().position) +
                std::abs((int)(lhs.get_inserted_sequence().size() - rhs.get_inserted_sequence().size())));
    }
    else
    {
        return std::numeric_limits<int>::infinity();
    }        
}

std::vector<Cluster> hierarchical_clustering_method(std::vector<Junction> const & junctions,
                                                    double clustering_cutoff)
{
    auto partitions = partition_junctions(junctions);
    std::vector<Cluster> clusters{};
    for (std::vector<Junction> partition : partitions)
    {
        size_t partition_size = partition.size();
        if (partition_size < 2)
        {
            clusters.emplace_back(partition);
            continue;
        }
        // Compute condensed distance matrix (upper triangle of the full distance matrix)
        double* distmat = new double[(partition_size * (partition_size - 1)) / 2];
        int k, i, j;
        for (i = k = 0; i < partition_size; ++i) {
            for (j = i + 1; j< partition_size; ++j) {
                // Compute distance between junctions i and j  
                distmat[k] = junction_distance(partition[i], partition[j]);
                ++k;
            }
        }

        // Perform hierarchical clustering
        // `height` is filled with cluster distance for each step
        // `merge` contains dendrogram
        int* merge = new int[2 * (partition_size - 1)];
        double* height = new double[partition_size - 1];
        hclust_fast(partition_size, distmat, HCLUST_METHOD_AVERAGE, merge, height);

        // Fill labels[i] with cluster label of junction i.
        // Clustering is stopped at step with cluster distance >= clustering_cutoff
        int* labels = new int[partition_size];
        cutree_cdist(partition_size, merge, height, clustering_cutoff, labels);

        std::unordered_map<int, std::vector<Junction>> label_to_junctions {};
        for (int i = 0; i < partition_size; ++i)
        {
            if (label_to_junctions.find(labels[i]) != label_to_junctions.end())
            {
                label_to_junctions[labels[i]].push_back(partition[i]);
            }
            else{
                label_to_junctions.emplace(labels[i], std::vector{partition[i]});
            }
        }

        // Add new clusters
        for (auto & [lab, jun] : label_to_junctions )
        {
            std::sort(jun.begin(), jun.end());
            clusters.emplace_back(jun);
        }

        // Free memory
        delete[] distmat;
        delete[] merge;
        delete[] height;
        delete[] labels;
    }
    std::sort(clusters.begin(), clusters.end());
    return clusters;
}
