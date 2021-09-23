#include <iostream>

#include <seqan3/alphabet/cigar/cigar.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>

#include "iGenVar.hpp"
#include "variant_detection/bam_functions.hpp"          // check the sam flags
#include "variant_detection/snp_indel_detection.hpp"    // detect_snp_and_indel

void detect_snp_and_indel(std::filesystem::path const & reads_filename)
{
    // Get the header information and set the necessary fields.
    using sam_fields = seqan3::fields<seqan3::field::flag,       // 2: FLAG
                                      seqan3::field::ref_id,     // 3: RNAME
                                      seqan3::field::ref_offset, // 4: POS
                                      seqan3::field::mapq,       // 5: MAPQ
                                      seqan3::field::cigar>;     // 6: CIGAR

    seqan3::sam_file_input reads_file{reads_filename, sam_fields{}};
    auto const & header = reads_file.header();

    // Check that the file is sorted.
    if (header.sorting != "coordinate")
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};

    // Prepare one activity vector per reference sequence.
    std::vector<std::vector<unsigned>> activity{};
    activity.resize(header.ref_id_info.size()); // allocate the number of reference sequences
    for (size_t idx = 0; idx < activity.size(); ++idx)
        activity[idx].resize(std::get<0>(header.ref_id_info[idx])); // allocate the length of the reference sequence

    for (auto && record : reads_file)
    {
        seqan3::sam_flag const flag = record.flag();                            // 2: FLAG
        int32_t const ref_id        = record.reference_id().value_or(-1);       // 3: RNAME
        int32_t ref_pos             = record.reference_position().value_or(-1); // 4: POS

        // Skip reads with certain properties.
        if (hasFlagUnmapped(flag) ||
            hasFlagSecondary(flag) ||
            hasFlagDuplicate(flag) ||
            record.mapping_quality() < 20 ||
            ref_id < 0 ||
            ref_pos < 0)
        {
            continue;
        }

        // Track the current position within the read.
        unsigned read_pos = 0;

        // Step through the CIGAR string.
        for (auto && [length, operation] : record.cigar_sequence())
        {
            // Skip long variants.
            if (length >= 30) // TODO: use min_length parameter
            {
                read_pos += length;
                ref_pos += length;
                continue;
            }

            // Case distinction for cigar elements.
            seqan3::cigar::operation const op = operation;
            switch (op.to_char())
            {
                case 'I':
                case 'S': // insertion or soft clip
                {
                    activity[ref_id][ref_pos] += length;
                    read_pos += length;
                    break;
                }
                case 'D': // deletion
                {
                    for (unsigned idx = 0; idx < length; ++idx)
                        ++activity[ref_id][ref_pos + idx];
                    ref_pos += length;
                    break;
                }
                default: // match or mismatch
                { // TODO: ATTENTION, mismatches are currently ignored and should still be verified with the reference
                    ref_pos += length;
                    read_pos += length;
                }
            }
        }
    }

    // Print the raw activity profiles.
    for (auto const & av : activity)
        seqan3::debug_stream << av << std::endl;
}
