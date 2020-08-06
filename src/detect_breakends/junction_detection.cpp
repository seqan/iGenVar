#include "detect_breakends/junction_detection.hpp"

using seqan3::operator""_cigar_op;
using seqan3::operator""_tag;

// ToDo: unused:
// std::unordered_map<std::string, int32_t> construct_ref_id_map(const std::deque<std::string> & ref_ids)
// {
//     std::unordered_map<std::string, int32_t> ref_id_map{};
//     int32_t index = 0;
//     for (std::string ref_name : ref_ids)
//     {
//         ref_id_map.insert({ref_name, index});
//         index++;
//     }
//     return ref_id_map;
// }

template <class Container>
void split_string(const std::string& str, Container& cont, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}

void retrieve_aligned_segments(std::string sa_string, std::vector<aligned_segment> & aligned_segments)
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
            aligned_segments.push_back(aligned_segment{ref_name, pos, orientation, cigar_vector, mapq});
        }
    }
}

void analyze_aligned_segments(const std::vector<aligned_segment> & aligned_segments,
                              std::vector<junction> & junctions,
                              std::string & read_name)
{
    for(size_t i = 1; i<aligned_segments.size(); i++)
    {
        aligned_segment current = aligned_segments[i-1];
        aligned_segment next = aligned_segments[i];
        int32_t distance_on_read = next.get_query_start() - current.get_query_end();
        // Neither gap nor overlap on read
        if (distance_on_read >= -10 && distance_on_read <= 10)
        {
            breakend mate1{current.ref_name,
                            current.orientation == strand::forward ? current.get_reference_end()
                                                                   : current.get_reference_start(),
                            current.orientation,
                            sequence_type::reference};
            breakend mate2{next.ref_name,
                            next.orientation == strand::forward ? next.get_reference_start()
                                                                : next.get_reference_end(),
                            next.orientation,
                            sequence_type::reference};
            junction new_junction{std::move(mate1), std::move(mate2), read_name};
            seqan3::debug_stream << "BND: " << new_junction << "\n";
            junctions.push_back(std::move(new_junction));
        }
    }
}

void analyze_cigar(std::vector<seqan3::cigar> & cigar_string,
                   std::vector<junction> & junctions,
                   std::vector<seqan3::dna5_vector> & insertions,
                   std::string chromosome,
                   int32_t query_start_pos,
                   seqan3::dna5_vector & query_sequence,
                   int32_t min_length,
                   seqan3::sequence_file_output<> & insertion_file,
                   std::string & read_name)
{
    // Step through CIGAR string and store current position in reference and read
    int32_t pos_ref = query_start_pos;
    int32_t pos_read = 0;

    // Stores the index of the current read in the insertion allele output file (or -1 if current read has not been added yet)
    int32_t insertion_allele_id {-1};

    for (seqan3::cigar & pair : cigar_string)
    {
        using seqan3::get;
        int32_t length = get<0>(pair);
        seqan3::cigar_op operation = get<1>(pair);
        if (operation == 'M'_cigar_op || operation == '='_cigar_op || operation == 'X'_cigar_op)
        {
            pos_ref += length;
            pos_read += length;
        }
        else if (operation == 'I'_cigar_op)
        {
            if (length > min_length)
            {
                if (insertion_allele_id < 0)
                {
                    insertion_allele_id = insertions.size();
                    std::string insertion_allele_name{"allele_" + std::to_string(insertion_allele_id)};
                    insertion_file.emplace_back(query_sequence, insertion_allele_name);
                    insertions.push_back(query_sequence);
                }
                // Insertions cause two junctions ( (1) from the reference to the read and (2) back from the read to the reference )
                junction new_junction1{breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                       breakend{std::to_string(insertion_allele_id),
                                                pos_read,
                                                strand::forward,
                                                sequence_type::read},
                                       read_name};
                seqan3::debug_stream << "INS1: " << new_junction1 << "\n";
                junctions.push_back(std::move(new_junction1));
                junction new_junction2{breakend{std::to_string(insertion_allele_id),
                                                pos_read + length,
                                                strand::forward,
                                                sequence_type::read},
                                       breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                       read_name};
                seqan3::debug_stream << "INS2: " << new_junction2 << "\n";
                junctions.push_back(std::move(new_junction2));
            }
            pos_read += length;
        }
        else if (operation == 'D'_cigar_op)
        {
            if (length > min_length)
            {
                // Deletions cause one junction from its start to its end
                junction new_junction{breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                      breakend{chromosome, pos_ref + length, strand::forward, sequence_type::reference},
                                      read_name};
                seqan3::debug_stream << "DEL: " << new_junction << "\n";
                junctions.push_back(std::move(new_junction));
            }
            pos_ref += length;
        }
        else if (operation == 'S'_cigar_op)
        {
            pos_read += length;
        }
    }
}

void detect_junctions_in_alignment_file(const std::filesystem::path & alignment_file_path,
                                        const std::filesystem::path & insertion_file_path)
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

    seqan3::alignment_file_input alignment_file{alignment_file_path, my_fields{}};
    // Open output file for insertion alleles
    seqan3::sequence_file_output insertion_file{insertion_file_path};

    // Store junctions, insertion_alleles and number of good alignments
    std::vector<junction> junctions{};
    std::vector<seqan3::dna5_vector> insertion_alleles{};
    uint16_t num_good = 0;

    for (auto & rec : alignment_file)
    {
        std::string query_name = seqan3::get<seqan3::field::id>(rec);
        int32_t ref_id = seqan3::get<seqan3::field::ref_id>(rec).value_or(0);
        int32_t pos = seqan3::get<seqan3::field::ref_offset>(rec).value_or(0);
        seqan3::sam_flag const flag = seqan3::get<seqan3::field::flag>(rec);    // uint16_t enum
        uint8_t const mapq = seqan3::get<seqan3::field::mapq>(rec);
        auto cigar = seqan3::get<seqan3::field::cigar>(rec);
        auto seq = seqan3::get<seqan3::field::seq>(rec);
        auto tags = seqan3::get<seqan3::field::tags>(rec);
        auto header_ptr = seqan3::get<seqan3::field::header_ptr>(rec);
        auto ref_ids = header_ptr->ref_ids();
        std::string ref_name = ref_ids[ref_id];

        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20)
        {
            // seqan3::debug_stream << "Skipped flag " << flag << std::endl;
        }
        else
        {
            // Detect junctions from CIGAR string
            analyze_cigar(cigar, junctions, insertion_alleles, ref_name, pos, seq, 30, insertion_file, query_name);
            // Detect junctions from SA tag (primary alignments only)
            if (!hasFlagSupplementary(flag))
            {
                std::string sa_tag = tags.get<"SA"_tag>();
                if (!sa_tag.empty())
                {
                    std::vector<aligned_segment> aligned_segments{};
                    auto strand = (hasFlagReverseComplement(flag) ? strand::reverse : strand::forward);
                    aligned_segments.push_back(aligned_segment{ref_name,  pos, strand, cigar, mapq});
                    retrieve_aligned_segments(sa_tag, aligned_segments);
                    std::sort(aligned_segments.begin(), aligned_segments.end());
                    analyze_aligned_segments(aligned_segments, junctions, query_name);
                }
            }

            num_good++;
            if (num_good % 1000 == 0)
            {
                seqan3::debug_stream << num_good << " good alignments" << std::endl;
            }
        }
    }
    std::sort(junctions.begin(), junctions.end());
    for (junction elem : junctions)
    {
        std::cout << elem << '\n';
    }
    seqan3::debug_stream << "Done. Found " << junctions.size() << " junctions." << '\n';
}
