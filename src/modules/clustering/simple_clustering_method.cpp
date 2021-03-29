#include "modules/clustering/simple_clustering_method.hpp"  // for the simple clustering method

#include <seqan3/core/debug_stream.hpp>

void simple_clustering_method(std::vector<Junction> const & junctions, std::vector<Cluster> & clusters)
{
    if (junctions.size() > 0)
    {
        std::vector<Junction> current_cluster_members = {junctions[0]};
        size_t i = 1;
        while (i < junctions.size())
        {
            if (junctions[i] == current_cluster_members.back())
            {
                current_cluster_members.push_back(junctions[i]);
            }
            else
            {
                clusters.emplace_back(std::move(current_cluster_members));
                current_cluster_members = {junctions[i]};
            }
            ++i;
        }
        clusters.emplace_back(std::move(current_cluster_members));
    }
    else
    {
        seqan3::debug_stream << "No junctions found...\n";
    }
}
