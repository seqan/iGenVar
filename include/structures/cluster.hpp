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
    //!\}

    Cluster(std::vector<Junction> members) : members{std::move(members)}
    {

    }
    size_t get_cluster_size() const;

    Breakend get_average_mate1() const;

    Breakend get_average_mate2() const;
};

template <typename stream_t>
inline stream_t operator<<(stream_t && stream, Cluster const & clust)
{
    stream << clust.get_average_mate1() << '\t'
           << clust.get_average_mate2() << '\t'
           << clust.get_cluster_size();
    return stream;
}
