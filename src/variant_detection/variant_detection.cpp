#include "variant_detection/variant_detection.hpp"

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>         // SAM/BAM support (seqan3::sam_file_input)

#include "modules/sv_detection_methods/analyze_cigar_method.hpp"    // for the split read method
#include "modules/sv_detection_methods/analyze_sa_tag_method.hpp"   // for the cigar string method
#include "variant_detection/bam_functions.hpp"                      // for hasFlag* functions

using seqan3::operator""_tag;

void detect_junctions_in_long_reads_sam_file(std::vector<Junction> & junctions,
                                             std::filesystem::path const & alignment_long_reads_file_path,
                                             std::vector<detection_methods> const & methods,
                                             clustering_methods const & clustering_method,
                                             refinement_methods const & refinement_method,
                                             uint64_t const min_var_length)
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
    seqan3::sam_file_input alignment_long_reads_file{alignment_long_reads_file_path, my_fields{}};

    // Check that the file is sorted before proceeding.
    if (alignment_long_reads_file.header().sorting != "coordinate")
    {
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};
    }
    uint16_t num_good = 0;

    for (auto & record : alignment_long_reads_file)
    {
        std::string const query_name        = record.id();                              // 1: QNAME
        seqan3::sam_flag const flag         = record.flag();                            // 2: FLAG
        int32_t const ref_id                = record.reference_id().value_or(-1);       // 3: RNAME
        int32_t const ref_pos               = record.reference_position().value_or(-1); // 4: POS
        uint8_t const mapq                  = record.mapping_quality();                 // 5: MAPQ
        std::vector<seqan3::cigar> cigar    = record.cigar_sequence();                  // 6: CIGAR
        auto const seq                      = record.sequence();                        // 10:SEQ
        auto tags                           = record.tags();
        auto const header_ptr               = record.header_ptr();
        auto const ref_ids = header_ptr->ref_ids();

        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20 ||
            ref_id < 0 || ref_pos < 0)
            continue;

        std::string const ref_name = ref_ids[ref_id];
        for (detection_methods method : methods) {
            switch (method)
            {
                case detection_methods::cigar_string: // Detect junctions from CIGAR string
                    analyze_cigar(query_name,
                                  ref_name,
                                  ref_pos,
                                  cigar,
                                  seq,
                                  junctions,
                                  min_var_length);
                    break;
                case detection_methods::split_read:     // Detect junctions from split read evidence (SA tag,
                    if (!hasFlagSupplementary(flag))    //                                  primary alignments only)
                    {
                        std::string const sa_tag = tags.get<"SA"_tag>();
                        if (!sa_tag.empty())
                        {
                            analyze_sa_tag(query_name, flag, ref_name, ref_pos, mapq, cigar, seq, sa_tag, junctions);
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
    std::sort(junctions.begin(), junctions.end());
}
