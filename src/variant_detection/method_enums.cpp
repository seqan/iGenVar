#include "variant_detection/method_enums.hpp"

std::unordered_map<std::string, detection_methods> enumeration_names(detection_methods)
{
return std::unordered_map<std::string, detection_methods>{{"0", detection_methods::cigar_string},
                                                          {"cigar_string", detection_methods::cigar_string},
                                                          {"1", detection_methods::split_read},
                                                          {"split_read", detection_methods::split_read},
                                                          {"2", detection_methods::read_pairs},
                                                          {"read_pairs", detection_methods::read_pairs},
                                                          {"3", detection_methods::read_depth},
                                                          {"read_depth", detection_methods::read_depth}};
};

std::unordered_map<std::string, clustering_methods> enumeration_names(clustering_methods)
{
return std::unordered_map<std::string,
                            clustering_methods>{{"0", clustering_methods::simple_clustering},
                                                {"simple_clustering",
                                                clustering_methods::simple_clustering},
                                                {"1", clustering_methods::hierarchical_clustering},
                                                {"hierarchical_clustering",
                                                clustering_methods::hierarchical_clustering},
                                                {"2", clustering_methods::self_balancing_binary_tree},
                                                {"self_balancing_binary_tree",
                                                clustering_methods::self_balancing_binary_tree},
                                                {"3", clustering_methods::candidate_selection_based_on_voting},
                                                {"candidate_selection_based_on_voting",
                                                clustering_methods::candidate_selection_based_on_voting}};
};

std::unordered_map<std::string, refinement_methods> enumeration_names(refinement_methods)
{
return std::unordered_map<std::string, refinement_methods>{{"0", refinement_methods::no_refinement},
                                                           {"no_refinement",
                                                           refinement_methods::no_refinement},
                                                           {"1", refinement_methods::sViper_refinement_method},
                                                           {"sViper_refinement_method",
                                                           refinement_methods::sViper_refinement_method},
                                                           {"2", refinement_methods::sVirl_refinement_method},
                                                           {"sVirl_refinement_method",
                                                           refinement_methods::sVirl_refinement_method}};
};
