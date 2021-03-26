#include "structures/cluster.hpp"

#include <cmath>        // for std::round
#include <stdexcept>    // for std::runtime_error

size_t Cluster::get_cluster_size() const
{
    return members.size();
}

Breakend Cluster::get_average_mate1() const
{
    std::string seq_name{};
    uint64_t sum_positions = 0;
    strand orientation{};
    // Iterate through members of the cluster
    for (uint32_t i = 0; i < members.size(); ++i)
    {
        Breakend mate1 = members[i].get_mate1();
        if (i == 0)
        {
            seq_name = mate1.seq_name;
            orientation = mate1.orientation;
        }
        else
        {
            // Make sure that all members of the cluster have matching sequence names and orientations
            if (mate1.seq_name != seq_name ||
                mate1.orientation != orientation)
            {
                throw std::runtime_error("Junctions with incompatible breakends were clustered together (different seq_name or orientation).");
            }
        }
        // Add up breakend positions aross all members
        sum_positions += mate1.position;
    }
    int32_t average_position = std::round(static_cast<double>(sum_positions) / members.size());
    Breakend average_breakend{seq_name, average_position, orientation};
    return average_breakend;
}

Breakend Cluster::get_average_mate2() const
{
    std::string seq_name{};
    uint64_t sum_positions = 0;
    strand orientation{};
    // Iterate through members of the cluster
    for (uint32_t i = 0; i < members.size(); ++i)
    {
        Breakend mate2 = members[i].get_mate2();
        if (i == 0)
        {
            seq_name = mate2.seq_name;
            orientation = mate2.orientation;
        }
        else
        {
            // Make sure that all members of the cluster have matching sequence types, names and orientations
            if (mate2.seq_name != seq_name ||
                mate2.orientation != orientation)
            {
                throw std::runtime_error("Junctions with incompatible breakends were clustered together (different seq_name or orientation).");
            }
        }
        // Add up breakend positions aross all members
        sum_positions += mate2.position;
    }
    int32_t average_position = std::round(static_cast<double>(sum_positions) / members.size());
    Breakend average_breakend{seq_name, average_position, orientation};
    return average_breakend;
}

int32_t Cluster::get_average_inserted_sequence_size() const
{
    uint64_t sum_sizes = 0;
    // Iterate through members of the cluster
    for (size_t i = 0; i < members.size(); ++i)
    {
        sum_sizes += members[i].get_inserted_sequence().size();
    }
    int32_t average_size = std::round(static_cast<double>(sum_sizes) / members.size());
    return average_size;
}