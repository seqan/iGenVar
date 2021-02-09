#include <seqan3/core/debug_stream.hpp>

#include "detect_breakends/bam_functions.hpp"   // for hasFlag* functions
#include "structures/aligned_segment.hpp"       // for struct AlignedSegment
#include "structures/junction.hpp"              // for class Junction

template <class Container>
void split_string(const std::string& str, Container& cont, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}

void retrieve_aligned_segments(std::string sa_string, std::vector<AlignedSegment> & aligned_segments)
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
            int32_t pos = std::stoi(fields[1]);
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
            std::tuple<std::vector<seqan3::cigar>, int32_t, int32_t> parsed_cigar = parse_cigar(cigar_field);
            std::vector<seqan3::cigar> cigar_vector = std::get<0>(parsed_cigar);
            int32_t mapq = std::stoi(fields[4]);
            aligned_segments.push_back(AlignedSegment{orientation, ref_name, pos, mapq, cigar_vector});
        }
    }
}

void analyze_aligned_segments(const std::vector<AlignedSegment> & aligned_segments,
                              std::vector<Junction> & junctions,
                              const std::string & read_name)
{
    for(size_t i = 1; i<aligned_segments.size(); i++)
    {
        AlignedSegment current = aligned_segments[i-1];
        AlignedSegment next = aligned_segments[i];
        int32_t distance_on_read = next.get_query_start() - current.get_query_end();
        // Neither gap nor overlap on read
        if (distance_on_read >= -10 && distance_on_read <= 10)
        {
            Breakend mate1{current.ref_name,
                            current.orientation == strand::forward ? current.get_reference_end()
                                                                   : current.get_reference_start(),
                            current.orientation,
                            sequence_type::reference};
            Breakend mate2{next.ref_name,
                            next.orientation == strand::forward ? next.get_reference_start()
                                                                : next.get_reference_end(),
                            next.orientation,
                            sequence_type::reference};
            Junction new_junction{std::move(mate1), std::move(mate2), read_name};
            seqan3::debug_stream << "BND: " << new_junction << "\n";
            junctions.push_back(std::move(new_junction));
        }
    }
}
