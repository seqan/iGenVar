#pragma once

//!\brief An enum for the different clustering methods.
enum detecting_methods
{
    cigar_string = 0,
    split_read = 1,
    read_pairs = 2,
    read_depth = 3,

    // Also add new methods to the default values in the argument parsers

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
