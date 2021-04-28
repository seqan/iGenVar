#include "variant_detection/variant_detection.hpp"

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>         // SAM/BAM support (seqan3::sam_file_input)

#include "modules/sv_detection_methods/analyze_cigar_method.hpp"    // for the split read method
#include "modules/sv_detection_methods/analyze_read_pair_method.hpp"// for the read pair method
#include "modules/sv_detection_methods/analyze_sa_tag_method.hpp"   // for the cigar string method
#include "variant_detection/bam_functions.hpp"                      // for hasFlag* functions

using seqan3::operator""_tag;

void detect_junctions_in_short_reads_sam_file(std::vector<Junction> & junctions,
                                              cmd_arguments const & args)
{
    // Open input alignment file
    using my_fields = seqan3::fields<seqan3::field::flag,       // 2: FLAG
                                     seqan3::field::ref_id,     // 3: RNAME
                                     seqan3::field::ref_offset, // 4: POS
                                     seqan3::field::mapq,       // 5: MAPQ
                                     seqan3::field::header_ptr>;

    seqan3::sam_file_input alignment_short_reads_file{args.alignment_short_reads_file_path, my_fields{}};

    // Check that the file is sorted before proceeding.
    if (alignment_short_reads_file.header().sorting != "coordinate")
    {
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};
    }
    uint16_t num_good = 0;

    for (auto & record : alignment_short_reads_file)
    {
        seqan3::sam_flag const flag         = record.flag();                            // 2: FLAG
        int32_t const ref_id                = record.reference_id().value_or(-1);       // 3: RNAME
        int32_t const ref_pos               = record.reference_position().value_or(-1); // 4: POS
        uint8_t const mapq                  = record.mapping_quality();                 // 5: MAPQ
        auto const header_ptr               = record.header_ptr();
        auto const ref_ids = header_ptr->ref_ids();

        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20 ||
            ref_id < 0 || ref_pos < 0)
            continue;

        std::string const ref_name = ref_ids[ref_id];
        for (detection_methods method : args.methods) {
            switch (method)
            {
                case detection_methods::cigar_string: // Detect junctions from CIGAR string
                    seqan3::debug_stream << "The cigar string method for short reads is not yet implemented.\n";
                    break;
                case detection_methods::split_read:     // Detect junctions from split read evidence (SA tag,
                    seqan3::debug_stream << "The split read method for short reads is not yet implemented.\n";
                    break;
                case detection_methods::read_pairs: // Detect junctions from read pair evidence
                    if (hasFlagMultiple(flag))
                    {
                        analyze_read_pair();
                    }
                    break;
                case detection_methods::read_depth: // Detect junctions from read depth evidence
                    seqan3::debug_stream << "The read depth method for short reads is not yet implemented.\n";
                    break;
            }
        }

        num_good++;
        if (num_good % 1000 == 0)
        {
            seqan3::debug_stream << num_good << " good alignments from short read file." << std::endl;
        }
    }
}

void detect_junctions_in_long_reads_sam_file(std::vector<Junction> & junctions,
                                             cmd_arguments const & args)
{
    // Open input alignment file
    using my_fields = seqan3::fields<seqan3::field::id,         // 1: QNAME
                                     seqan3::field::flag,       // 2: FLAG
                                     seqan3::field::ref_id,     // 3: RNAME
                                     seqan3::field::ref_offset, // 4: POS
                                     seqan3::field::mapq,       // 5: MAPQ
                                     seqan3::field::cigar,      // 6: CIGAR
                                     seqan3::field::seq,        // 10:SEQ
                                     seqan3::field::tags,
                                     seqan3::field::header_ptr>;
    seqan3::sam_file_input alignment_long_reads_file{args.alignment_long_reads_file_path, my_fields{}};

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
        seqan3::dna5_vector const seq       = record.sequence();                        // 10:SEQ
        auto tags                           = record.tags();
        auto const header_ptr               = record.header_ptr();
        auto const ref_ids = header_ptr->ref_ids();

        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20 ||
            ref_id < 0 || ref_pos < 0)
            continue;

        std::string const ref_name = ref_ids[ref_id];
        for (detection_methods method : args.methods) {
            switch (method)
            {
                case detection_methods::cigar_string: // Detect junctions from CIGAR string
                    analyze_cigar(query_name,
                                  ref_name,
                                  ref_pos,
                                  cigar,
                                  seq,
                                  junctions,
                                  args.min_var_length);
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
                case detection_methods::read_pairs: // There are no read pairs in long reads.
                    break;
                case detection_methods::read_depth: // Detect junctions from read depth evidence
                    seqan3::debug_stream << "The read depth method for long reads is not yet implemented.\n";
                    break;
            }
        }

        num_good++;
        if (num_good % 1000 == 0)
        {
            seqan3::debug_stream << num_good << " good alignments from long read file." << std::endl;
        }
    }
}
