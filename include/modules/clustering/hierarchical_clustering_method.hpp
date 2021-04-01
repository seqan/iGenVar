#pragma once

#include "structures/cluster.hpp"   // for class Cluster


/*! \brief Partition junctions by their distance on the reference genome.
 *         The junctions in each of the returned partitions are sorted even though
 *         the partitions themselves are not returned in a particular order.
 *
 * \param junctions - a vector of junctions (needs to be sorted)
 */
std::vector<std::vector<Junction>> partition_junctions(std::vector<Junction> const & junctions);

/*! \brief Sub-partition an existing partition based on the second mate of each junction.
 *         The junctions in each of the returned sub-partitions are sorted even though
 *         the sub-partitions themselves are not returned in a particular order.
 *
 * \param partition - a partition (i.e. a vector) of junctions (needs to be sorted)
 */
std::vector<std::vector<Junction>> split_partition_based_on_mate2(std::vector<Junction> const & partition);

/*! \brief Cluster junctions by an hierarchical clustering method.
 *
 * \param junctions - a sorted vector of junctions
 * \param clusters - a vector of clusters
 */
void hierarchical_clustering_method(std::vector<Junction> const & junctions, std::vector<Cluster> & clusters);
