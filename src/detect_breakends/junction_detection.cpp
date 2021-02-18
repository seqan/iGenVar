#include "detect_breakends/junction_detection.hpp"

using seqan3::operator""_cigar_op;
using seqan3::operator""_tag;

/*! \brief Splits a string by a given delimiter and stores substrings in a given container.
 *
 * \param str   string to split
 * \param cont  container for the splitted substrings
 * \param delim delimiter
 */
template <class Container>
void split_string(const std::string& str, Container& cont, char delim = ' ')
{
    std::stringstream ss(str);
    std::string token;
    while (std::getline(ss, token, delim)) {
        cont.push_back(token);
    }
}

/*! \brief Parse the SA tag from the SAM/BAM alignment of a chimeric/split-aligned read. Build
 *         [aligned_segments](\ref AlignedSegment), one for each alignment segment of the read.
 *
 * \param sa_string         "SA" tag string
 * \param aligned_segments  vector of [aligned_segments](\ref AlignedSegment).
 *
 * \details The SA tag describes the alignments of a chimeric read and is like a small SAM within a SAM file:
 *          "SA:Z:(rname,pos,strand,CIGAR,mapQ,NM;)+"
 *          Each element (in parentheses) represents one alignment segment of the chimeric alignment formatted as
 *          a colon-delimited list.
 *          We add all segments to our candidate list `aligned_segments` and examine them in the following function
 *          `analyze_aligned_segments()`.
 *          For more information about this tag, see the
 *          ([Map Optional Fields Specification](https://github.com/samtools/hts-specs/blob/master/SAMtags.pdf)).
 */
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

/*! \brief Build junctions out of aligned_segments.
 *
 * \param aligned_segments  vector of [aligned_segments](\ref AlignedSegment).
 * \param junctions         vector for storing junctions
 * \param read_name         QNAME field of the SAM/BAM file
 */
void analyze_aligned_segments(const std::vector<AlignedSegment> & aligned_segments,
                              std::vector<Junction> & junctions,
                              std::string & read_name)
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

/*! \brief This function steps through the CIGAR string and stores junctions with their position in reference and read.
 *
 * \param chromosome        RNAME field of the SAM/BAM file
 * \param read_name         QNAME field of the SAM/BAM file
 * \param query_start_pos   POS field of the SAM/BAM file
 * \param cigar_string      CIGAR field of the SAM/BAM file
 * \param query_sequence    SEQ field of the SAM/BAM file
 * \param junctions         vector for storing junctions
 * \param insertions        vector for storing insertion_alleles
 * \param min_length        minimum length of variants to detect (default 30 bp)
 * \param insertion_file    output file for insertion alleles
 *
 * \details This function steps through the CIGAR string and stores junctions with their position in reference and read.
 *          We distinguish 4 cases of CIGAR operation characters:
 *          1. M, =, X: For matches and mismatches (Alignment column containing two letters. This could contain two
 *                      different letters (mismatch) or two identical letters), we step through ref and read.
 *          2. I:       For insertons (gap in the reference sequence) -> Insertions cause two junctions ( (1) from the
 *                      reference to the read and (2) back from the read to the reference ).
 *          3. D:       For deletions (gap in the query sequence) -> Deletions cause one junction from its start to its
 *                      end.
 *          4. S:       For soft clipped letters we step through the read. These are segments of the query sequence that
 *                      do not appear in the alignment. The full-length query sequence is given in the SEQ field of the
 *                      SAM record.
 *          Other CIGAR operations: H, N, P are skipped (H: hard clipping sequences are not present in the SEQ, N:
 *          skipped region representing an intron, P: padding consumes neither the query nor the reference).
 *          The junctions found are stored in the given `junctions` vector.
 *          For more information see the
 *          ([Map Format Specification](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf)) page 8.
 */
void analyze_cigar(std::string chromosome,
                   std::string & read_name,
                   int32_t query_start_pos,
                   std::vector<seqan3::cigar> & cigar_string,
                   seqan3::dna5_vector & query_sequence,
                   std::vector<Junction> & junctions,
                   std::vector<seqan3::dna5_vector> & insertions,
                   int32_t min_length,
                   seqan3::sequence_file_output<> & insertion_file)
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
        else if (operation == 'I'_cigar_op) // I: Insertion (gap in the reference sequence).
        {
            if (length >= min_length)
            {
                if (insertion_allele_id < 0)
                {
                    insertion_allele_id = insertions.size();
                    std::string insertion_allele_name{"allele_" + std::to_string(insertion_allele_id)};
                    insertion_file.emplace_back(query_sequence, insertion_allele_name);
                    insertions.push_back(query_sequence);
                }
                // Insertions cause two junctions ( (1) from the reference to the read and (2) back from the read to the reference )
                Junction new_junction1{Breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                       Breakend{std::to_string(insertion_allele_id),
                                                pos_read,
                                                strand::forward,
                                                sequence_type::read},
                                       read_name};
                seqan3::debug_stream << "INS1: " << new_junction1 << "\n";
                junctions.push_back(std::move(new_junction1));
                Junction new_junction2{Breakend{std::to_string(insertion_allele_id),
                                                pos_read + length,
                                                strand::forward,
                                                sequence_type::read},
                                       Breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                       read_name};
                seqan3::debug_stream << "INS2: " << new_junction2 << "\n";
                junctions.push_back(std::move(new_junction2));
            }
            pos_read += length;
        }
        else if (operation == 'D'_cigar_op)
        {
            if (length >= min_length)
            {
                // Deletions cause one junction from its start to its end
                Junction new_junction{Breakend{chromosome, pos_ref, strand::forward, sequence_type::reference},
                                      Breakend{chromosome, pos_ref + length, strand::forward, sequence_type::reference},
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
        else // other possible cigar operations: H, N, P
        {
            // seqan3::debug_stream << "Unhandled operation " << operation << std::endl;
        }

    }
}

/*! \brief Detects junctions between distant genomic positions by analyzing an alignment file (sam/bam). The detected
 *         junctions are printed on stdout and insertion alleles are stored in a fasta file.
 * \cond
 * \param alignment_file_path input file - path to the sam/bam file
 * \param insertion_file_path output file - path for the fasta file
 * \param methods - list of methods for detecting junctions (1: cigar_string,
 *                                                           2: split_read,
 *                                                           3: read_pairs,
 *                                                           4: read_depth)
 * \param clustering_method method for clustering junctions (0: simple_clustering
 *                                                           1: hierarchical_clustering,
 *                                                           2: self-balancing_binary_tree,
 *                                                           3: candidate_selection_based_on_voting)
 * \param refinement_method method for refining breakends (0: no_refinement,
 *                                                         1: sViper_refinement_method,
 *                                                         2: sVirl_refinement_method)
 * \param min_var_length - minimum length of variants to detect (default 30 bp)
 * \endcond
 *
 * \details Detects junctions from the CIGAR strings and supplementary alignment tags of read alignment records.
 *          In the iteration over all reads, we first sort out unmapped alignments, secondary alignments, duplicates
 *          and alignments with low mapping quality. You can find information about all flags at the
 *          ([Map Format Specification](https://github.com/samtools/hts-specs/blob/master/SAMv1.pdf)) page 7.
 *          Then, the CIGAR string of all remaining alignments is analyzed.
 *          For primary alignments, the SA tag is analyzed additionally yielding information on supplementary alignments
 *          of a split read.
 *          More details on this in the associated function `retrieve_aligned_segments()`.
 */
void detect_junctions_in_alignment_file(const std::filesystem::path & alignment_file_path,
                                        const std::filesystem::path & insertion_file_path,
                                        const std::vector<uint8_t> methods,
                                        const clustering_methods clustering_method,
                                        const refinement_methods refinement_method,
                                        const uint64_t min_var_length)
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
    std::vector<Junction> junctions{};
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
            for (uint8_t method : methods) {
                switch (method)
                {
                    case 1: // Detect junctions from CIGAR string
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
                    case 2: // Detect junctions from split read evidence (SA tag, primary alignments only)
                        if (!hasFlagSupplementary(flag))
                        {
                            std::string sa_tag = tags.get<"SA"_tag>();
                            if (!sa_tag.empty())
                            {
                                std::vector<AlignedSegment> aligned_segments{};
                                auto strand = (hasFlagReverseComplement(flag) ? strand::reverse : strand::forward);
                                aligned_segments.push_back(AlignedSegment{strand, ref_name, pos, mapq, cigar});
                                retrieve_aligned_segments(sa_tag, aligned_segments);
                                std::sort(aligned_segments.begin(), aligned_segments.end());
                                analyze_aligned_segments(aligned_segments, junctions, query_name);
                            }
                        }
                        break;
                    case 3: // Detect junctions from read pair evidence
                        seqan3::debug_stream << "The read pair method is not yet implemented.\n";
                        break;
                        // continue;
                    case 4: // Detect junctions from read depth evidence
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
    }
    std::sort(junctions.begin(), junctions.end());

    seqan3::debug_stream << "Start clustering...\n";

    std::vector<Cluster> clusters{};
    switch (clustering_method)
    {
        case 0: // simple_clustering
            {
                std::vector<Junction> current_cluster_members = {junctions[0]};
                int i = 1;
                while (i < junctions.size())
                {
                    if (junctions[i] == current_cluster_members.back())
                    {
                        current_cluster_members.push_back(junctions[i]);
                    }
                    else
                    {
                        clusters.emplace_back(current_cluster_members);
                        current_cluster_members = {junctions[i]};
                    }
                    ++i;
                }
                clusters.emplace_back(current_cluster_members);
            }
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

    for (Cluster const & elem : clusters)
    {
        std::cout << elem << '\n';
    }
}
