#include "variant_detection/variant_detection.hpp"

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>         // SAM/BAM support (seqan3::sam_file_input)

#include "modules/sv_detection_methods/analyze_cigar_method.hpp"    // for the split read method
#include "modules/sv_detection_methods/analyze_read_pair_method.hpp"// for the read pair method
#include "modules/sv_detection_methods/analyze_sa_tag_method.hpp"   // for the cigar string method
#include "variant_detection/bam_functions.hpp"                      // for hasFlag* functions

using seqan3::operator""_tag;

std::deque<std::string> read_header_information(auto & alignment_file,
                                                std::map<std::string, int32_t> & references_lengths)
{
    // Check that the file is sorted before proceeding.
    if (alignment_file.header().sorting != "coordinate")
    {
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};
    }

    // Get the information from \@SQ tag, more precise the values of the SN and LN tags
    std::deque<std::string> const ref_ids = alignment_file.header().ref_ids();
    std::vector<std::tuple<int32_t, std::string>> ref_id_info = alignment_file.header().ref_id_info;

    size_t i = 0;
    for(std::string const & ref_id : ref_ids)
    {
        int32_t ref_length = std::get<0>(ref_id_info[i]);
        if (references_lengths.find(ref_id) != references_lengths.end())
        {
            if (references_lengths[ref_id] != ref_length)
            {
                std::cerr << "Warning: The reference id " << ref_id << " was found twice in the input files with "
                          << "different length: " << references_lengths[ref_id] << " and " << ref_length << '\n';
            }
        } else {
            references_lengths.emplace(ref_id, ref_length);
        }
        ++i;
    }

    return ref_ids;
}

void detect_junctions_in_short_reads_sam_file(std::vector<Junction> & junctions,
                                              std::map<std::string, int32_t> & references_lengths,
                                              cmd_arguments const & args)
{
    // Open input alignment file
    using my_fields = seqan3::fields<seqan3::field::flag,       // 2: FLAG
                                     seqan3::field::ref_id,     // 3: RNAME
                                     seqan3::field::ref_offset, // 4: POS
                                     seqan3::field::mapq>;      // 5: MAPQ

    // Set number of decompression threads
    seqan3::contrib::bgzf_thread_count = args.threads;
    seqan3::sam_file_input alignment_short_reads_file{args.alignment_short_reads_file_path, my_fields{}};

    std::deque<std::string> const ref_ids = read_header_information(alignment_short_reads_file, references_lengths);
    uint32_t num_good = 0;

    for (auto & record : alignment_short_reads_file)
    {
        seqan3::sam_flag const flag         = record.flag();                            // 2: FLAG
        int32_t const ref_id                = record.reference_id().value_or(-1);       // 3: RNAME
        int32_t const ref_pos               = record.reference_position().value_or(-1); // 4: POS
        uint8_t const mapq                  = record.mapping_quality();                 // 5: MAPQ
        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20 ||
            ref_id < 0 || ref_pos < 0)
            continue;

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
                                             std::map<std::string, int32_t> & references_lengths,
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
                                     seqan3::field::tags>;

    // Set number of decompression threads
    seqan3::contrib::bgzf_thread_count = args.threads;
    seqan3::sam_file_input alignment_long_reads_file{args.alignment_long_reads_file_path, my_fields{}};

    std::deque<std::string> const ref_ids = read_header_information(alignment_long_reads_file, references_lengths);
    uint32_t num_good = 0;

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
                            analyze_sa_tag(query_name,
                                           flag,
                                           ref_name,
                                           ref_pos,
                                           mapq,
                                           cigar,
                                           seq,
                                           sa_tag,
                                           args,
                                           junctions);
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
