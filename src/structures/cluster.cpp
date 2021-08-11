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

size_t Cluster::get_common_tandem_dup_count() const
{
    size_t sum_counts = 0;
    size_t amount_zero_counts = 0;
    // Iterate through members of the cluster
    for (size_t i = 0; i < members.size(); ++i)
    {
        size_t current_count = members[i].get_tandem_dup_count();
        if (current_count == 0)
            ++amount_zero_counts;
        else
            sum_counts += members[i].get_tandem_dup_count();
    }
    // TODO (irallia 12.08.21): This 3 is free choosen and other values could be tested.
    // If two thirds of the junctions have a 0 tandem_dup_count, than its probably no tandem duplication.
    if (amount_zero_counts > std::round(members.size() / 3.0))
        return 0;
    else
        return std::round(static_cast<double>(sum_counts) / members.size());
}

size_t Cluster::get_average_inserted_sequence_size() const
{
    size_t sum_sizes = 0;
    // Iterate through members of the cluster
    for (size_t i = 0; i < members.size(); ++i)
    {
        sum_sizes += members[i].get_inserted_sequence().size();
    }
    size_t average_size = std::round(static_cast<double>(sum_sizes) / members.size());
    return average_size;
}

std::vector<Junction> Cluster::get_members() const
{
    return members;
}

bool operator<(Cluster const & lhs, Cluster const & rhs)
{
    return lhs.get_average_mate1() != rhs.get_average_mate1()
           ? lhs.get_average_mate1() < rhs.get_average_mate1()
           : lhs.get_average_mate2() != rhs.get_average_mate2()
              ? lhs.get_average_mate2() < rhs.get_average_mate2()
              : lhs.get_common_tandem_dup_count() != rhs.get_common_tandem_dup_count()
                ? lhs.get_common_tandem_dup_count() < rhs.get_common_tandem_dup_count()
                : lhs.get_average_inserted_sequence_size() < rhs.get_average_inserted_sequence_size();
}

bool operator==(Cluster const & lhs, Cluster const & rhs)
{
    std::vector<Junction> lhs_members = lhs.get_members();
    std::sort(lhs_members.begin(), lhs_members.end());
    std::vector<Junction> rhs_members = rhs.get_members();
    std::sort(rhs_members.begin(), rhs_members.end());
    return (lhs_members == rhs_members);
}
