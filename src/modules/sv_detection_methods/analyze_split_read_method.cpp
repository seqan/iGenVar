#include "iGenVar.hpp"              // for global variable gVerbose
#include "modules/sv_detection_methods/analyze_split_read_method.hpp"

#include <seqan3/core/debug_stream.hpp>

using seqan3::operator""_dna5;

template <class Container>
void split_string(std::string const & str, Container & cont, char const delim)
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}

void retrieve_aligned_segments(std::string const & sa_string, std::vector<AlignedSegment> & aligned_segments)
{
    std::vector<std::string> sa_tags{};
    split_string(sa_string, sa_tags, ';');
    for (std::string sa_tag : sa_tags)
    {
        std::vector<std::string> fields {};
        split_string(sa_tag, fields, ',');
        if (fields.size() == 6)
        {
            std::string ref_name = fields[0];
            // Decrement by 1 because position in SA tag is stored as string and 1-based unlike other coordinates
            int32_t pos = std::stoi(fields[1]) - 1;
            strand orientation;
            if (fields[2] == "+")
            {
                orientation = strand::forward;
            }
            else if (fields[2] == "-")
            {
                orientation = strand::reverse;
            }
            else
            {
                continue;
            }
            std::string cigar_field = fields[3];
            std::vector<seqan3::cigar> cigar_vector = std::get<0>(seqan3::detail::parse_cigar(cigar_field));
            int32_t mapq = std::stoi(fields[4]);
            aligned_segments.push_back(AlignedSegment{orientation, ref_name, pos, mapq, cigar_vector});
        }
        else
        {
            seqan3::debug_stream << "Your SA tag has a wrong format (wrong amount of parameters): " << sa_tag << '\n';
        }
    }
}

std::pair<int32_t, int32_t> get_mate_positions(AlignedSegment current, AlignedSegment next, int32_t distance_on_read)
{
    int32_t mate1_pos;
    if (current.orientation == strand::forward)
        mate1_pos = current.get_reference_end();
    else
        mate1_pos = current.get_reference_start();
    int32_t mate2_pos;
    if (next.orientation == strand::forward)
    {
        // Correct position of mate 2 for overlapping alignment segments:
        // Trim alignment of `next` segment at the start to remove overlap
        if (distance_on_read < 0)
            mate2_pos = next.get_reference_start() - distance_on_read;
        else
            mate2_pos = next.get_reference_start();
    }
    else
    {
        // Correct position of mate 2 for overlapping alignment segments:
        // Trim alignment of `next` segment at the end to remove overlap
        if (distance_on_read < 0)
            mate2_pos = next.get_reference_end() + distance_on_read;
        else
            mate2_pos = next.get_reference_end();
    }
    return std::pair(mate1_pos, mate2_pos);
}

void analyze_aligned_segments(std::vector<AlignedSegment> const & aligned_segments,
                              std::vector<Junction> & junctions,
                              seqan3::dna5_vector const & query_sequence,
                              std::string const & read_name,
                              int32_t const min_length,
                              int32_t const max_overlap)
{
    bool last_element_was_tandem_duplication = false;
    for (size_t i = 1; i < aligned_segments.size(); i++)
    {
        AlignedSegment current = aligned_segments[i-1];
        AlignedSegment next = aligned_segments[i];
        int32_t distance_on_read = next.get_query_start() - current.get_query_end();
        // Check that the overlap between two consecutive alignment segments
        // of the read is lower than the given threshold
        if (distance_on_read >= -max_overlap)
        {
            auto [ mate1_pos, mate2_pos ] = get_mate_positions(current, next, distance_on_read);
            int32_t distance_on_ref = mate2_pos - mate1_pos - 1;
            // Check that the two consecutive alignment segments either
            // map to different reference sequences (e.g. translocation, interspersed duplication),
            // have a different orientation (e.g. inversion)
            // have a large distance on the reference (e.g. deletion, inversion, tandem duplication), or
            // have a large distance on the read (e.g. insertion)
            if (current.ref_name != next.ref_name ||        // 1: Segments not on same chromosome -> translocation
                std::abs(distance_on_ref) >= min_length ||  // 3: Distance on reference >= minimum SV size -> 4
                distance_on_read >= min_length)
            {
                Breakend mate1{current.ref_name, mate1_pos, current.orientation};
                Breakend mate2{next.ref_name, mate2_pos, next.orientation};
                size_t tandem_dup_count = 0;
                // if novel inserted sequence
                if (current.ref_name == next.ref_name && distance_on_read > 0)
                {
                    last_element_was_tandem_duplication = false;
                    auto inserted_bases = query_sequence | seqan3::views::slice(current.get_query_end(),
                                                                                next.get_query_start());
                    junctions.emplace_back(mate1, mate2, inserted_bases, tandem_dup_count, read_name);
                    if (gVerbose)
                        seqan3::debug_stream << "INS: " << junctions.back() << "\n";
                }
                else
                {
                    // If overlapping aligned segments
                    //                         |-DUP:TANDEM-|
                    // ref ----------------------------------------------
                    //                 ||||||||||||||||||||||
                    // current_segment ----------------------
                    //                         ||||||||||||||||||||||
                    // next_segment            ----------------------
                    if (current.ref_name == next.ref_name &&
                        current.get_reference_start() <= next.get_reference_start() &&
                        next.get_reference_start() < current.get_reference_end())
                    {
                        // If multiple overlapping aligned segments (next pair is also duplicated)
                        //                         |-DUP:TANDEM-|
                        // ref ----------------------------------------------
                        //                 ||||||||||||||||||||||
                        // current_segment ----------------------
                        //                         ||||||||||||||
                        // next_segment            --------------
                        //                         ||||||||||||||||||||||
                        // next_but_one_segment    ----------------------
                        if(last_element_was_tandem_duplication)
                        {
                            tandem_dup_count = junctions.back().get_tandem_dup_count() + 1;
                            // Replace last element
                            junctions.back() = Junction{mate2, mate1, ""_dna5, tandem_dup_count, read_name};
                        }
                        else // first pair of a tandem duplication
                        {
                            tandem_dup_count+=2;
                            junctions.emplace_back(mate2, mate1, ""_dna5, tandem_dup_count, read_name);
                        }
                        last_element_was_tandem_duplication = true;
                        if (gVerbose)
                            seqan3::debug_stream << "DUP:TANDEM: " << junctions.back() << "\n";
                    }
                    else
                    { // Else TRA or DUP or INV
                        last_element_was_tandem_duplication = false;
                        junctions.emplace_back(mate1, mate2, ""_dna5, tandem_dup_count, read_name);
                        if (gVerbose)
                            seqan3::debug_stream << "BND: " << junctions.back() << "\n";
                    }

                }
            }
        }
    }
}

void analyze_sa_tag(std::string const & query_name,
                    seqan3::sam_flag const & flag,
                    std::string const & ref_name,
                    int32_t const pos,
                    uint8_t const mapq,
                    std::vector<seqan3::cigar> const & cigar,
                    seqan3::dna5_vector const & seq,
                    std::string const & sa_tag,
                    cmd_arguments const & args,
                    std::vector<Junction> & junctions)
{
    std::vector<AlignedSegment> aligned_segments{};
    strand strand = (hasFlagReverseComplement(flag) ? strand::reverse : strand::forward);
    aligned_segments.push_back(AlignedSegment{strand, ref_name, pos, mapq, cigar});
    retrieve_aligned_segments(sa_tag, aligned_segments);
    // sort by query start, query end, mapping quality (in this order):
    std::sort(aligned_segments.begin(), aligned_segments.end());
    analyze_aligned_segments(aligned_segments, junctions, seq, query_name, args.min_var_length, args.max_overlap);
}
