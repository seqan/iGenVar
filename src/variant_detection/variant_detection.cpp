#include "variant_detection/variant_detection.hpp"

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>         // SAM/BAM support
#include <seqan3/io/sequence_file/output.hpp>   // FASTA support

#include "misc/statistics.hpp"                                      // for calculating average read depth
#include "modules/clustering/simple_clustering_method.hpp"          // for the simple clustering method
#include "modules/sv_detection_methods/analyze_cigar_method.hpp"    // for the split read method
#include "modules/sv_detection_methods/analyze_sa_tag_method.hpp"   // for the cigar string method
#include "structures/cluster.hpp"                                   // for class Cluster
#include "variant_detection/bam_functions.hpp"                      // for hasFlag* functions
#include "variant_detection/variant_output.hpp"

using seqan3::operator""_tag;


void detect_variants_in_alignment_file(const std::filesystem::path & alignment_file_path,
                                       const std::filesystem::path & insertion_file_path,
                                       const std::vector<detection_methods> & methods,
                                       const clustering_methods & clustering_method,
                                       const refinement_methods & refinement_method,
                                       const uint64_t & min_var_length,
                                       const std::filesystem::path & output_file_path,
                                       const uint64_t & sample_size)
{
    // Open input alignment file
    using my_fields = seqan3::fields<seqan3::field::id,
                                     seqan3::field::ref_id,
                                     seqan3::field::ref_offset,
                                     seqan3::field::flag,
                                     seqan3::field::mapq,
                                     seqan3::field::cigar,
                                     seqan3::field::seq,
                                     seqan3::field::tags,
                                     seqan3::field::header_ptr>;

    seqan3::sam_file_input alignment_file{alignment_file_path, my_fields{}};
    // Open output file for insertion alleles
    seqan3::sequence_file_output insertion_file{insertion_file_path};

    // Store junctions, insertion_alleles and number of good alignments
    std::vector<Junction> junctions{};
    std::vector<seqan3::dna5_vector> insertion_alleles{};
    uint16_t num_good = 0;

    // Obtain list of random chromosomes and positions to get read depth from.
    std::forward_list<std::pair<int32_t, int32_t>> rand_sample = get_random_positions(sample_size, alignment_file.header());
    uint32_t avg_depth{0};
    uint32_t cur_depth{0};
    uint32_t num_pos{0};

    for (auto & rec : alignment_file)
    {
        const std::string query_name            = seqan3::get<seqan3::field::id>(rec);                      // 1: QNAME
        const seqan3::sam_flag flag             = seqan3::get<seqan3::field::flag>(rec);                    // 2: FLAG
        const int32_t ref_id                    = seqan3::get<seqan3::field::ref_id>(rec).value_or(-1);     // 3: RNAME
        const int32_t pos                       = seqan3::get<seqan3::field::ref_offset>(rec).value_or(-1);  // 4: POS
        const uint8_t mapq                      = seqan3::get<seqan3::field::mapq>(rec);                    // 5: MAPQ
        std::vector<seqan3::cigar> cigar        = seqan3::get<seqan3::field::cigar>(rec);                   // 6: CIGAR
        const auto seq                          = seqan3::get<seqan3::field::seq>(rec);                     // 10:SEQ
        auto tags                               = seqan3::get<seqan3::field::tags>(rec);
        const auto header_ptr                   = seqan3::get<seqan3::field::header_ptr>(rec);
        const auto ref_ids = header_ptr->ref_ids();

        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20 ||
            ref_id < 0 || pos < 0)
            continue;

        const std::string ref_name = ref_ids[ref_id];
        uint32_t read_length{get_read_length(cigar)};
        // If there are positions to be sampled, then sample them.
        if (!rand_sample.empty() && sample_size > 0)
        {
            if (ref_id == std::get<0>(rand_sample.front()) &&           // Correct chromosome.
                pos <= std::get<1>(rand_sample.front()) &&              // Read starts before or at the position.
                std::get<1>(rand_sample.front()) <= pos + read_length)  // Read ends after or at position.
            {
                ++cur_depth;
            }
            else if (ref_id > std::get<0>(rand_sample.front()) ||       // Record is at next chromosome
                     (ref_id == std::get<0>(rand_sample.front()) &&     // Record is at same chromosome BUT
                      pos > std::get<1>(rand_sample.front())))          // in a later position.
            {
                ++num_pos;
                avg_depth += cur_depth;
                cur_depth = 0;
                rand_sample.pop_front();
            }
        }
        for (uint8_t method : methods) {
            switch (method)
            {
                case detection_methods::cigar_string: // Detect junctions from CIGAR string
                    analyze_cigar(query_name,
                                  ref_name,
                                  pos,
                                  cigar,
                                  seq,
                                  junctions,
                                  insertion_alleles,
                                  min_var_length,
                                  insertion_file);
                    break;
                case detection_methods::split_read:     // Detect junctions from split read evidence (SA tag,
                    if (!hasFlagSupplementary(flag))    //                                  primary alignments only)
                    {
                        const std::string sa_tag = tags.get<"SA"_tag>();
                        if (!sa_tag.empty())
                        {
                            analyze_sa_tag(query_name, flag, ref_name, pos, mapq, cigar, sa_tag, junctions);
                        }
                    }
                    break;
                case detection_methods::read_pairs: // Detect junctions from read pair evidence
                    seqan3::debug_stream << "The read pair method is not yet implemented.\n";
                    break;
                    // continue;
                case detection_methods::read_depth: // Detect junctions from read depth evidence
                    seqan3::debug_stream << "The read depth method is not yet implemented.\n";
                    break;
                    // continue;
            }
        }

        num_good++;
        if (num_good % 1000 == 0)
        {
            seqan3::debug_stream << num_good << " good alignments" << std::endl;
        }
    }
    if (sample_size == 0 || num_pos == 0) // User requested no sampling, or all of the random positions generated were not covered.
    {
        avg_depth = 0;
    }
    else
    {
        avg_depth = avg_depth/num_pos;
    }
    std::sort(junctions.begin(), junctions.end());

    seqan3::debug_stream << "Start clustering...\n";

    std::vector<Cluster> clusters{};
    switch (clustering_method)
    {
        case 0: // simple_clustering
            simple_clustering_method(junctions, clusters);
            break;
        case 1: // hierarchical clustering
            seqan3::debug_stream << "The hierarchical clustering method is not yet implemented\n";
            break;
        case 2: // self-balancing_binary_tree,
            seqan3::debug_stream << "The self-balancing binary tree clustering method is not yet implemented\n";
            break;
        case 3: // candidate_selection_based_on_voting
            seqan3::debug_stream << "The candidate selection based on voting clustering method is not yet implemented\n";
            break;
    }

    seqan3::debug_stream << "Done with clustering. Found " << clusters.size() << " junction clusters.\n";

    switch (refinement_method)
    {
        case 0: // no refinement
            seqan3::debug_stream << "No refinement was selected.\n";
            break;
        case 1: // sViper_refinement_method
            seqan3::debug_stream << "The sViper refinement method is not yet implemented\n";
            break;
        case 2: // sVirl_refinement_method
            seqan3::debug_stream << "The sVirl refinement method is not yet implemented\n";
            break;
    }

    find_and_output_variant(clusters, output_file_path);
}
