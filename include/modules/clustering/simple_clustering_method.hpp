#pragma once

#include "structures/cluster.hpp"   // for class Cluster

/*! \brief This method clusters junctions by merging identical junctions into a single cluster object.
 *
 * \param[in]       junctions   - a vector of junctions
 * \param[in, out]  clusters    - a vector of clusters
 */
void simple_clustering_method(std::vector<Junction> & junctions, std::vector<Cluster> & clusters);
