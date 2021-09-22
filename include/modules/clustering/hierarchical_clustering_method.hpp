#pragma once

#include "structures/cluster.hpp"   // for class Cluster

/*! \brief Partition junctions by their distance on the reference genome.
 *         The returned partitions contain junctions meeting the following criteria:
 *         a) all junctions in a partition connect the same reference sequences,
 *         b) all junctions in a partition have the same orientations, and
 *         c) the distance between corresponding mates of two neighboring junctions is
 *            smaller or equal than the given partition_max_distance.
 *         The junctions in each of the returned partitions are sorted even though
 *         the partitions themselves are not returned in a particular order.
 *
 * \param[in] junctions - a vector of junctions (needs to be sorted)
 * \param[in] partition_max_distance - maximum distance between junctions in the same partition
 *
 * \returns Returns unordered partitions of sorted junctions.
 */
std::vector<std::vector<Junction>> partition_junctions(std::vector<Junction> const & junctions,
                                                       int32_t const partition_max_distance);

/*! \brief Sub-partition an existing partition based on the second mate of each junction.
 *         The junctions in each of the returned sub-partitions are sorted even though
 *         the sub-partitions themselves are not returned in a particular order.
 *
 * \param[in] partition - a partition (i.e. a vector) of junctions (needs to be sorted
 *                        by their second mates)
 * \param[in] partition_max_distance - maximum distance between junctions in the same partition
 *
 * \returns Returns unordered sub-partitions of sorted junctions based on the second mate of each junction.
 */
std::vector<std::vector<Junction>> split_partition_based_on_mate2(std::vector<Junction> const & partition,
                                                                  int32_t const partition_max_distance);

/*! \brief Compute the distance between two junctions.
 *         For two junctions that connect the same reference sequences and have the same
 *         orientations, the distance is the sum of a) the distance between the first mates,
 *         b) the distance between the second mates, and c) the absolute size difference
 *         of the inserted sequences. For other junctions, the distance takes the maximal value.
 *
 * \param[in] lhs - left side junction
 * \param[in] rhs - right side junction
 *
 * \returns Returns distance between two junctions.
 */
double junction_distance(Junction const & lhs, Junction const & rhs);

/*! \brief Cluster junctions by an hierarchical clustering method.
 *         The returned clusters and the junctions in each returned cluster are sorted.
 *
 * \param[in] junctions - a vector of junctions (needs to be sorted)
 * \param[in] partition_max_distance - maximum distance between junctions in the same partition
 * \param[in] clustering_cutoff - distance cutoff for clustering
 *
 * \returns Returns sorted clusters of sorted junctions.
 *
 * \details For the algorithms we use the library hclust.
 * \see https://lionel.kr.hs-niederrhein.de/~dalitz/data/hclust/ (last access 01.06.2021).
 */
std::vector<Cluster> hierarchical_clustering_method(std::vector<Junction> const & junctions,
                                                    int32_t const partition_max_distance,
                                                    double clustering_cutoff);
