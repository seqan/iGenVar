#include "modules/clustering/hierarchical_clustering_method.hpp"  // for the hierarchical clustering method

#include <seqan3/core/debug_stream.hpp>

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

void hierarchical_clustering_method(std::vector<Junction> const & junctions, std::vector<Cluster> & clusters)
{
    auto partitions = partition_junctions(junctions);
    for (std::vector<Junction> partition : partitions)
    {
        //TODO (eldariont): replace with hierarchical clustering algorithm 
        //                  that clusters junctions in each partition
        //                  e.g. from https://github.com/cdalitz/hclust-cpp
        clusters.push_back(std::move(partition));
    }
}
