#pragma once

#include <vector>

#include "structures/junction.hpp"  // for class Junction

class Cluster
{
private:
    std::vector<Junction> members{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr Cluster()                    = default; //!< Defaulted.
    Cluster(Cluster const &)               = default; //!< Defaulted.
    Cluster(Cluster &&)                    = default; //!< Defaulted.
    Cluster & operator=(Cluster const &)   = default; //!< Defaulted.
    Cluster & operator=(Cluster &&)        = default; //!< Defaulted.
    ~Cluster()                             = default; //!< Defaulted.

    Cluster(std::vector<Junction> members) : members{std::move(members)}
    {

    }
    //!\}

    // \brief Returns the number of members in the cluster.
    size_t get_cluster_size() const;

    /*! \brief Returns the average first mate of all cluster members.
    *          All cluster members are required to have identical sequence names and orientations for their first mate.
    *          To produce the average, the average first mate's position of all cluster members is computed.
    */
    Breakend get_average_mate1() const;

    /*! \brief Returns the average second mate of all cluster members.
    *          All cluster members are required to have identical sequence names and orientations for their second mate.
    *          To produce the average, the average second mate's position of all cluster members is computed.
    */
    Breakend get_average_mate2() const;

    // \brief Returns the average length of the inserted sequences of all cluster members.
    int32_t get_average_inserted_sequence_size() const;
};

template <typename stream_t>
inline stream_t operator<<(stream_t && stream, Cluster const & clust)
{
    stream << clust.get_average_mate1() << '\t'
           << clust.get_average_mate2() << '\t'
           << clust.get_cluster_size() << '\t'
           << clust.get_average_inserted_sequence_size();
    return stream;
}
