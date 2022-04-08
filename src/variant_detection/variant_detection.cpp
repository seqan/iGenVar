#include "variant_detection/variant_detection.hpp"

#include <fcntl.h>
#include <unistd.h>

#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>         // SAM/BAM support (seqan3::sam_file_input)

#include "modules/sv_detection_methods/analyze_cigar_method.hpp"        // for the split read method
#include "modules/sv_detection_methods/analyze_read_pair_method.hpp"    // for the read pair method
#include "modules/sv_detection_methods/analyze_split_read_method.hpp"   // for the cigar string method
#include "variant_detection/bam_functions.hpp"                          // for hasFlag* functions

#include "cereal/types/memory.hpp"
#include "cereal/types/vector.hpp"

using seqan3::operator""_tag;

// SAM fields for input file
using my_fields = seqan3::fields<seqan3::field::id,         // 1: QNAME
                                 seqan3::field::flag,       // 2: FLAG
                                 seqan3::field::ref_id,     // 3: RNAME
                                 seqan3::field::ref_offset, // 4: POS
                                 seqan3::field::mapq,       // 5: MAPQ
                                 seqan3::field::cigar,      // 6: CIGAR
                                 seqan3::field::seq,        // 10:SEQ
                                 seqan3::field::tags>;

std::deque<std::string> read_header_information(auto & alignment_file,
                                                std::map<std::string, int32_t> & references_lengths)
{
    // Check that the file is sorted before proceeding.
    if (alignment_file.header().sorting != "coordinate")
    {
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};
    }

    // Get the information from \@SQ tag, more precise the values of the SN and LN tags
    std::deque<std::string> ref_ids = alignment_file.header().ref_ids();
    std::vector<std::tuple<int32_t, std::string>> ref_id_info = alignment_file.header().ref_id_info;

    size_t i = 0;
    for (std::string & ref_id : ref_ids)
    {
        // As the chromosome names can differ between "1" and "chr1", we would distinguish same SVs from different
        // input files if the chromosome naming is different. Thus we add the "chr" prefix.
        if (!ref_id.starts_with("chr")) {
            // Add "chr" prefix
            ref_id = "chr" + ref_id;
        }
        int32_t ref_length = std::get<0>(ref_id_info[i]);
        if (references_lengths.find(ref_id) != references_lengths.end())
        {
            if (references_lengths[ref_id] != ref_length)
            {
                std::cerr << "Warning: The reference id " << ref_id << " was found twice in the input files with "
                          << "different length: " << references_lengths[ref_id] << " and " << ref_length << '\n';
            }
        }
        else
        {
            references_lengths.emplace(ref_id, ref_length);
        }
        ++i;
    }

    return ref_ids;
}

void safe_sync_rename(std::filesystem::path const & tmp_file_path, std::filesystem::path const & file_path) {
    int fd = open(tmp_file_path.string().c_str(), O_APPEND);
    fsync(fd);
    close(fd);
    std::filesystem::rename(tmp_file_path, file_path);
}

std::vector<std::unique_ptr<bamit::IntervalNode>> load_or_create_index(std::filesystem::path const & input_path)
{
    std::filesystem::path bamit_index_file_path{input_path};
    bamit_index_file_path += ".bit";
    std::vector<std::unique_ptr<bamit::IntervalNode>> node_list{};
    if (std::filesystem::exists(bamit_index_file_path))
    {
        std::ifstream in{bamit_index_file_path, std::ios_base::binary | std::ios_base::in};
        cereal::BinaryInputArchive ar(in);
        bamit::read(node_list, ar);
        in.close();
    }
    else
    {
        std::filesystem::path bamit_index_file_tmp_path{input_path};
        bamit_index_file_tmp_path += ".bit.tmp";
        std::ofstream out{bamit_index_file_tmp_path, std::ios_base::binary | std::ios_base::out};
        cereal::BinaryOutputArchive ar(out);
        seqan3::sam_file_input input_sam{input_path, my_fields{}};
        node_list = bamit::index(input_sam);
        bamit::write(node_list, ar);
        out.close();
        // If a run of iGenVar is aborted during BAMIT creation, an incorrect index is stored, which cannot be read when
        // iGenVar is called again. Therefore, we write the index to a tmp file and save it with the correct file name
        // when finished.
        safe_sync_rename(bamit_index_file_tmp_path, bamit_index_file_path);
    }
    return node_list;
}

void detect_junctions_in_short_reads_sam_file([[maybe_unused]] std::vector<Junction> & junctions,
                                              std::map<std::string, int32_t> & references_lengths,
                                              cmd_arguments const & args)
{
    // Open input alignment file
    seqan3::sam_file_input alignment_short_reads_file{args.alignment_short_reads_file_path, my_fields{}};

    std::deque<std::string> const ref_ids = read_header_information(alignment_short_reads_file, references_lengths);
    uint32_t num_good = 0;

    // Load bamit index, or create index if it doesn't exist.
    std::vector<std::unique_ptr<bamit::IntervalNode>> bamit_index = load_or_create_index(args.alignment_short_reads_file_path);

    for (auto & record : alignment_short_reads_file)
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

        std::string const & ref_name = ref_ids[ref_id];

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
        if (gVerbose)
        {
            ++num_good;
            if (num_good % 100000 == 0)
            {
                seqan3::debug_stream << num_good << " good alignments from short read file." << std::endl;
            }
        }
    }
}

void detect_junctions_in_long_reads_sam_file(std::vector<Junction> & junctions,
                                             std::map<std::string, int32_t> & references_lengths,
                                             cmd_arguments const & args)
{
    // Open input alignment file
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

        std::string const & ref_name = ref_ids[ref_id];

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
                case detection_methods::read_pairs:
                    // There are no read pairs in long reads.
                    break;
                case detection_methods::read_depth: // Detect junctions from read depth evidence
                    seqan3::debug_stream << "The read depth method for long reads is not yet implemented.\n";
                    break;
            }
        }

        if (gVerbose)
        {
            ++num_good;
            if (num_good % 100000 == 0)
            {
                seqan3::debug_stream << num_good << " good alignments from long read file." << std::endl;
            }
        }
    }
}
