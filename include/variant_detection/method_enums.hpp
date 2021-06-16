#pragma once

#include <string>           // for std::string
#include <unordered_map>    // for std::unordered_map

//!\brief An enum for the different clustering methods.
enum detection_methods
{
    cigar_string = 0,
    split_read = 1,
    read_pairs = 2,
    read_depth = 3,

    // Also add new methods to the default values in the argument parser

    //!\cond
    // ATTENTION: Must always be the last item; will be used to determine the number of ids.
    SIZE //!< Determines the size of the enum.
    //!\endcond
};

//!\brief An enum for the different clustering methods.
enum clustering_methods
{
    simple_clustering = 0,
    hierarchical_clustering = 1,
    self_balancing_binary_tree = 2,
    candidate_selection_based_on_voting = 3
};

//!\brief An enum for the different refinement methods.
enum refinement_methods
{
    no_refinement = 0,
    sViper_refinement_method = 1,
    sVirl_refinement_method = 2
};

/*! \brief Specialise a mapping from an identifying string to the respective value of your type detection_methods. With
 *         the help of this function, you're able to call ./detect_breackends with --method 0 and --method cigar_string
 *         and get the same result.
 */
std::unordered_map<std::string, detection_methods> enumeration_names(detection_methods);

/*! \brief Specialise a mapping from an identifying string to the respective value of your type clustering_methods. With
 *         the help of this function, you're able to call ./detect_breackends with --clustering_method 0 and
 *         --clustering_method simple_clustering and get the same result.
 */
std::unordered_map<std::string, clustering_methods> enumeration_names(clustering_methods);

/*! \brief Specialise a mapping from an identifying string to the respective value of your type refinement_methods. With
 *         the help of this function, you're able to call ./detect_breackends with --refinement_method 0 and
 *         --refinement_method no_refinement and get the same result.
 */
std::unordered_map<std::string, refinement_methods> enumeration_names(refinement_methods);
