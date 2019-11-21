#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include <seqan3/io/sequence_file/all.hpp>  // FASTA support
#include <seqan3/io/alignment_file/all.hpp> // SAM/BAM support

#include "bam_functions.hpp"
#include "junction.hpp"
#include "aligned_segment.hpp"

using namespace seqan3;

template <class Container>
void split_string(const std::string& str, Container& cont, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}

void retrieve_aligned_segments(std::string sa_string, std::vector<aligned_segment> & aligned_segments, std::unordered_map<std::string, int32_t> ref_id_map)
{
    std::vector<std::string> sa_tags{};
    split_string(sa_string, sa_tags, ';');
    for (std::string sa_tag : sa_tags)
    {
        std::vector<std::string> fields {};
        split_string(sa_tag, fields, ',');
        if (fields.size() == 6)
        {
            std::string rname = fields[0];
            int32_t ref_id = ref_id_map[rname];
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
            std::tuple<std::vector<cigar>, int32_t, int32_t> parsed_cigar = parse_cigar(cigar_field);
            std::vector<cigar> cigar_vector = std::get<0>(parsed_cigar);
            int32_t mapq = std::stoi(fields[4]);
            aligned_segments.push_back(aligned_segment{ref_id, pos, orientation, cigar_vector, mapq});
        }
    }
}

void analyze_aligned_segments(const std::vector<aligned_segment> & aligned_segments, std::vector<junction> & junctions)
{
    for(size_t i = 1; i<aligned_segments.size(); i++)
    {
        aligned_segment current = aligned_segments[i-1];
        aligned_segment next = aligned_segments[i];
        int32_t distance_on_read = next.get_query_start() - current.get_query_end();
        // Neither gap nor overlap on read
        if (distance_on_read >= -10 && distance_on_read <= 10)
        {
            breakend mate1{current.ref_id,
                            current.orientation == strand::forward ? current.get_reference_end()
                                                                   : current.get_reference_start(),
                            current.orientation,
                            sequence_type::reference};
            breakend mate2{next.ref_id,
                            next.orientation == strand::forward ? next.get_reference_start()
                                                                : next.get_reference_end(),
                            next.orientation,
                            sequence_type::reference};
            junction new_junction{std::move(mate1), std::move(mate2)};
            debug_stream << "BND: " << new_junction << "\n";
            junctions.push_back(std::move(new_junction));
        }
    }
}

void analyze_cigar(std::vector<cigar> & cigar_string, std::vector<junction> & junctions, std::vector<dna5_vector> & insertions, int32_t chromosome, int32_t query_start_pos, dna5_vector & query_sequence, int32_t min_length, sequence_file_output<> & insertion_file)
{
    // Step through CIGAR string and store current position in reference and read
    int32_t pos_ref = query_start_pos;
    int32_t pos_read = 0;

    // Stores the index of the current read in the insertion allele output file (or -1 if current read has not been added yet)
    int32_t insertion_allele_id {-1};

    for (cigar & pair : cigar_string)
    {
        int32_t length = get<0>(pair);
        cigar_op operation = get<1>(pair);
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
                junction new_junction1{breakend{chromosome, query_start_pos, strand::forward, sequence_type::reference},
                                       breakend{insertion_allele_id, pos_read, strand::forward, sequence_type::read}};
                debug_stream << "INS1: " << new_junction1 << "\n";
                junctions.push_back(std::move(new_junction1));
                junction new_junction2{breakend{insertion_allele_id, pos_read + length, strand::forward, sequence_type::read},
                                       breakend{chromosome, query_start_pos, strand::forward, sequence_type::reference}};
                debug_stream << "INS2: " << new_junction2 << "\n";
                junctions.push_back(std::move(new_junction2));
            }
            pos_read += length;
        }
        else if (operation == 'D'_cigar_op)
        {
            if (length > min_length)
            {
                // Deletions cause one junction from its start to its end
                junction new_junction{breakend{chromosome, query_start_pos, strand::forward, sequence_type::reference},
                                      breakend{chromosome, query_start_pos + length, strand::forward, sequence_type::reference}};
                debug_stream << "DEL: " << new_junction << "\n";
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

std::unordered_map<std::string, int32_t> construct_ref_id_map(const std::deque<std::string> & ref_ids)
{
    std::unordered_map<std::string, int32_t> ref_id_map{};
    int32_t index = 0;
    for (std::string ref_name : ref_ids)
    {
        ref_id_map.insert({ref_name, index});
        index++;
    }
    return ref_id_map;
}


void detect_junctions_in_alignment_file(const std::filesystem::path & alignment_file_path, const std::filesystem::path & insertion_file_path)
{
    // Open input alignment file
    using my_fields = fields<field::ID,
                             field::REF_ID,
                             field::REF_OFFSET,
                             field::FLAG,
                             field::MAPQ,
                             field::CIGAR,
                             field::SEQ,
                             field::TAGS,
                             field::HEADER_PTR>;

    alignment_file_input alignment_file{alignment_file_path, my_fields{}};

    // Open output file for insertion alleles
    sequence_file_output insertion_file{insertion_file_path};

    // Store junctions, insertion_alleles and number of good alignments
    std::vector<junction> junctions{};
    std::vector<dna5_vector> insertion_alleles{};
    uint16_t num_good = 0;
    bool ref_id_map_initialised = false;
    std::unordered_map<std::string, int32_t> ref_id_map;

    for (auto & rec : alignment_file)
    {
        std::string query_name = get<field::ID>(rec);
        int32_t ref_id = get<field::REF_ID>(rec).value_or(0);
        int32_t pos = get<field::REF_OFFSET>(rec).value_or(0);
        auto flag = get<field::FLAG>(rec);
        auto mapq = get<field::MAPQ>(rec);
        auto cigar = get<field::CIGAR>(rec);
        auto seq = get<field::SEQ>(rec);
        auto tags = get<field::TAGS>(rec);
        auto header_ptr = get<field::HEADER_PTR>(rec);
        auto ref_ids = header_ptr->ref_ids();
        std::string ref_name = ref_ids[ref_id];

        if (!ref_id_map_initialised)
        {
            ref_id_map = construct_ref_id_map(ref_ids);
            ref_id_map_initialised = true;
        }

        if (hasFlagUnmapped(flag) || hasFlagSecondary(flag) || hasFlagDuplicate(flag) || mapq < 20)
        {
            // debug_stream << "Skipped flag " << flag << std::endl;
        }
        else
        {
            // Detect junctions from CIGAR string
            analyze_cigar(cigar, junctions, insertion_alleles, ref_id, pos, seq, 30, insertion_file);
            // Detect junctions from SA tag (primary alignments only)
            if (!hasFlagSupplementary(flag))
            {
                std::string sa_tag = tags.get<"SA"_tag>();
                if (!sa_tag.empty())
                {
                    std::vector<aligned_segment> aligned_segments{};
                    aligned_segments.push_back(aligned_segment{ref_id, pos, hasFlagReverseComplement(flag) ? strand::reverse : strand::forward, cigar, mapq});
                    retrieve_aligned_segments(sa_tag, aligned_segments, ref_id_map);
                    std::sort(aligned_segments.begin(), aligned_segments.end());
                    analyze_aligned_segments(aligned_segments, junctions);
                }
            }

            num_good++;
            if (num_good % 1000 == 0)
            {
                debug_stream << num_good << " good alignments" << std::endl;
            }
        }
    }
    std::sort(junctions.begin(), junctions.end());
    for (junction elem : junctions)
    {
        std::cout << elem << '\n';
    }
    debug_stream << "Done. Found " << junctions.size() << " junctions." << '\n';
}

struct cmd_arguments
{
    std::filesystem::path alignment_file_path{};
    std::filesystem::path insertion_file_path{};
};

void initialize_argument_parser(argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "David Heller";
    parser.info.short_description = "Detect junctions in a read alignment file";
    parser.info.version = "0.0.1";
    parser.add_positional_option(args.alignment_file_path, "Input read alignments in SAM or BAM format.",
                                 input_file_validator{{"sam", "bam"}} );
    parser.add_positional_option(args.insertion_file_path, "Output file for insertion alleles",
                                 output_file_validator{{"fa", "fasta"}} );
}

int main(int argc, char ** argv)
{
    argument_parser myparser{"detectJunctions", argc, argv};        // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);
    try
    {
        myparser.parse();                                          // trigger command line parsing
    }
    catch (parser_invalid_argument const & ext)                     // catch user errors
    {
        debug_stream << "[Error] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
    detect_junctions_in_alignment_file(args.alignment_file_path, args.insertion_file_path);
    return 0;
}
