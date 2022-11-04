#include <iostream>
#include <numeric>
#include <sstream>

#include <seqan3/alignment/configuration/align_config_method.hpp>
#include <seqan3/alignment/configuration/align_config_output.hpp>
#include <seqan3/alignment/pairwise/align_pairwise.hpp>
#include <seqan3/alphabet/nucleotide/dna4.hpp>
#include <seqan3/core/debug_stream.hpp>
#include <seqan3/io/sam_file/input.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/utility/range/to.hpp>
#include <seqan3/utility/views/convert.hpp>
#include <seqan3/utility/views/elements.hpp>
#include <seqan3/utility/views/repeat.hpp>
#include <seqan3/utility/views/slice.hpp>
#include <seqan3/utility/views/zip.hpp>

#include <bamit/all.hpp>

#include "iGenVar.hpp"
#include "structures/debruijn_graph.hpp"
#include "variant_detection/bam_functions.hpp"          // check the sam flags
#include "variant_detection/snp_indel_detection.hpp"    // detect_snp_and_indel

Genome read_genome(std::filesystem::path const & genome_filename)
{
    Genome genome;
    seqan3::sequence_file_input sfi{genome_filename};
    for (auto & record : sfi)
    {
        genome.seqs.push_back(record.sequence());
        genome.names.push_back(record.id());
    }
    return genome;
}

/*!
 * \brief Extract activity from SAM records by counting indels and soft clips.
 * \param[in,out] activity       - The activity values for one reference genome.
 * \param[in]     cigar_sequence - The cigar string of the SAM record.
 * \param[in]     min_var_length - The length above which an indel/SNP is considered a variant.
 * \param[in]     ref_pos        - The start position of the alignment in the genome.
 */
void update_activity_for_record(std::vector<unsigned> & activity,
                                std::vector<seqan3::cigar> const & cigar_sequence,
                                uint64_t const min_var_length,
                                int32_t ref_pos)
{
    // Track the current position within the read.
    unsigned read_pos = 0;

    // Step through the CIGAR string.
    for (auto && [length, operation] : cigar_sequence)
    {
        // Skip long variants.
        if (length >= min_var_length)
        {
            read_pos += length;
            ref_pos += length;
            continue;
        }

        // Case distinction for cigar elements.
        switch (operation.to_char())
        {
            case 'I':
            case 'S': // insertion or soft clip
            {
                activity[ref_pos] += 1u;
                read_pos += length;
                break;
            }
            case 'D': // deletion
            {
                for (unsigned idx = 0; idx < length; ++idx)
                    ++activity[ref_pos + idx];
                ref_pos += length;
                break;
            }
            default: // ignore match, mismatch, hard clipping, padding, and skipped region
            { // TODO (joergi-w 30.09.2021) mismatches are ignored but should be verified with the reference
                ref_pos += length;
                read_pos += length;
            }
        }
    }
}

Activity analyze_activity(std::filesystem::path const & reads_filename, Genome const & genome, uint64_t min_var_length)
{
    Activity activity;
    activity.values.reserve(genome.seqs.size());
    activity.refmap.reserve(genome.seqs.size());

    // Get the header information and set the necessary fields.
    using sam_fields = seqan3::fields<seqan3::field::flag,       // 2: FLAG
                                      seqan3::field::ref_id,     // 3: RNAME
                                      seqan3::field::ref_offset, // 4: POS
                                      seqan3::field::mapq,       // 5: MAPQ
                                      seqan3::field::cigar>;     // 6: CIGAR

    seqan3::sam_file_input reads_file{reads_filename, sam_fields{}};
    auto & header = reads_file.header();

    // Check that the file is sorted.
    if (header.sorting != "coordinate")
        throw seqan3::format_error{"ERROR: Input file must be sorted by coordinate (e.g. samtools sort)"};

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

        // Check if we have already allocated the length of the reference sequence.
        if (activity.refmap.empty() || activity.refmap.back() != ref_id)
        {
            // Compensate genome sequences that do not appear as reference, if ref_id matches genome id.
            if (genome.seqs.size() == header.ref_ids().size())
            {
                while (activity.refmap.size() < static_cast<size_t>(ref_id))
                {
                    activity.refmap.push_back(-1); // -1 means unused
                    activity.values.emplace_back(0); // new vector of zero length
                }
            }

            // Treat the next genome sequence as current reference.
            activity.refmap.push_back(ref_id);
            activity.values.emplace_back(std::get<0>(header.ref_id_info[ref_id])); // new 0-vector of seq length

            // Some error handling.
            if (genome.seqs.size() < activity.refmap.size())
            {
                auto && rin = std::ranges::transform_view(activity.refmap,
                                                          [&header](int32_t idx){ return header.ref_ids()[idx]; });
                std::ostringstream oss;
                seqan3::debug_stream_type str{oss};
                str << "The genome file does not contain enough sequences: " << rin.size() << " expected: " << rin
                    << ", " << genome.names.size() << " found: " << genome.names;
                throw std::runtime_error(oss.str().c_str());
            }
            else if (activity.values.back().size() != genome.seqs[activity.refmap.size() - 1].size())
            {
                std::ostringstream oss;
                oss << "Mismatch of genome " << genome.names[activity.refmap.size() - 1] << " of length "
                    << genome.seqs[activity.refmap.size() - 1].size() << " and reference " << header.ref_ids()[ref_id]
                    << " of length " << activity.values.back().size() << ".";
                throw std::runtime_error(oss.str().c_str());
            }
        }

        // Fill the activity vector
        update_activity_for_record(activity.values.back(), record.cigar_sequence(), min_var_length, ref_pos);
    }
    return activity;
}

std::vector<ActiveRegions> get_active_regions(Activity activity, size_t window_width, size_t act_threshold)
{
    // TODO (joergi-w 31.05.2022) depend activity threshold on read coverage?
    // TODO (joergi-w 31.05.2022) another idea for window width: (min_var_length + 1) / 2;

    // Set up the result vector.
    std::vector<ActiveRegions> regions(activity.values.size());
    for (auto && [region, activ] : seqan3::views::zip(regions, activity.values))
    {
        size_t window_score = 0; // sum of the activities inside the window
        size_t pos;              // the current position in the genome sequence
        bool active = false;     // whether `pos` is inside an active region

        // Initialisation of a window of length min_var_length.
        for (pos = 0; pos < window_width && pos < activ.size(); ++pos)
            window_score += activ[pos];

        while (pos < activ.size())
        {
            // Check for start or end of an active region.
            if (!active && window_score >= act_threshold)
            {
                region.emplace_back(static_cast<int>(pos - window_width), 0);
                active = true;
            }
            else if (active && window_score < act_threshold)
            {
                std::get<1>(region.back()) = static_cast<int>(pos); // we report past-the-end position
                active = false;
            }

            // Advance the sliding window.
            window_score -= activ[pos - window_width];
            window_score += activ[pos];
            ++pos;
        }

        // Handle the case if the genome ends in an active region.
        if (active)
            std::get<1>(region.back()) = static_cast<int>(pos);
    }

    return regions;
}

std::vector<std::pair<float, seqan3::dna4_vector>> find_haplotypes(int32_t ref_id,
                                                                   seqan3::dna4_vector const & ref,
                                                                   std::pair<int, int> const & region,
                                                                   std::filesystem::path const & reads_filename)
{
    DeBruijnGraph graph;
    unsigned char len = 10; // kmer length is one of (10, 25, 32) -> SeqAn restricts len <= 32

    // try to create a graph as long as the graph is inviable, kmer length <= 32, and at least 1 kmer fits in the region
    while (!graph.is_viable() && len <= 32 && region.first + len - 1 < region.second - len + 1)
    {
        // Try to initialize the graph.
        if (graph.init_sequence(len, ref))
        {
            // start is the ref. position of the last character of the first kmer inside the region
            int32_t start = region.first + len - 1;
            // end is the ref. position of the first character of the last kmer inside the region
            int32_t end = region.second - len + 1;
            // from the loop condition we ensure that start < end

            // Open reads file to extract reads of regions.
            seqan3::sam_file_input reads_file{reads_filename};
            std::vector<std::unique_ptr<bamit::IntervalNode>> node_list = bamit::index(reads_file);

            // Get the file position of the first record matching the region.
            std::streamoff file_pos{-1};
            bamit::get_overlap_file_position(reads_file, node_list, {ref_id, start}, {ref_id, end}, file_pos);

            if (file_pos != -1)
            {   // Collect the reads as long as they start before the end of the region.
                auto w_fn = [end](auto & rec){ return rec.reference_position().value() < end; };
                // Keep only mapped reads. Throw away reads that end before the region starts.
                auto f_fn = [start](auto & rec){ return !bamit::unmapped(rec) && rec.reference_position().value() +
                                                 static_cast<int>(rec.sequence().size()) > start; };
                // Extract only the substring of the read that overlaps the region. Convert to dna4.
                auto t_fn = [region](auto & rec){ return rec.sequence() | seqan3::views::slice(std::max(0,
                                                  region.first - rec.reference_position().value()),
                                                  region.second - rec.reference_position().value())
                                                  | seqan3::views::convert<seqan3::dna4>; };
                auto && reads = reads_file
                                | std::views::take_while(w_fn)
                                | std::views::filter(f_fn)
                                | std::views::transform(t_fn);

                for (auto const & read : reads)
                    graph.add_read(seqan3::ranges::to<seqan3::dna4_vector>(read));
            }
        }
        len += len == 10 ? 15 : 7; // we try max 3 different kmer lengths (10, 25, 32)
    }
    if (graph.is_viable())
    {
        graph.prune(2);
        graph.export_dot_format(std::string{"dot/"} + std::to_string(region.first) + ".gv");
        auto haplotypes = graph.collect_haplotype_sequences();
        // Limit the number of haplotypes (they are sorted, so we choose the best ones).
        if (haplotypes.size() > 128)
            haplotypes.erase(haplotypes.begin() + 128, haplotypes.end());
        return haplotypes;
    }
    else
    {
        return std::vector<std::pair<float, seqan3::dna4_vector>>{};
    }
}

void store_snp(std::set<Junction> & junctions,
               seqan3::dna5_vector & bufR,
               seqan3::dna5_vector & bufH,
               int pos,
               std::string const & ref_name,
               float qual)
{
    if (bufR.empty() && bufH.empty())
        return;

    auto res = junctions.emplace(Breakend{ref_name, pos - 1 - static_cast<int>(bufR.size()), strand::forward},
                                 Breakend{ref_name, pos, strand::forward},
                                 bufH, bufR, 0, "", qual);
    if (!res.second)
    {
        qual += res.first->get_quality();
        junctions.erase(res.first);
        junctions.emplace(Breakend{ref_name, pos - 1 - static_cast<int>(bufR.size()), strand::forward},
                          Breakend{ref_name, pos, strand::forward},
                          bufH, bufR, 0, "", qual);
    }

//    DEBUG OUTPUT
//    if (bufR.empty()) // insertion
//    {
//        seqan3::debug_stream << "Found insertion at " << pos-1 << "-" << pos << ": " << bufH << " (" << qual << ")\n";
//    }
//    else if (bufH.empty()) // deletion
//    {
//        seqan3::debug_stream << "Found deletion at " << pos - 1 - bufR.size() << "-" << pos << ": " << bufR
//                             << " (" << qual << ")\n";
//    }
//    else // substitution
//    {
//        seqan3::debug_stream << "Found substitution at " << pos - 1 - bufR.size() << "-" << pos << ": "
//                             << bufR << " => " << bufH << " (" << qual << ")\n";
//    }

    bufR.clear();
    bufH.clear();
}

void call_snp(std::set<Junction> & junctions,
              seqan3::dna4_vector reference,
              std::vector<std::pair<float, seqan3::dna4_vector>> haplotypes,
              int region_start,
              std::string const & ref_name)
{
    auto config = seqan3::align_cfg::method_global{seqan3::align_cfg::free_end_gaps_sequence1_leading{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_leading{false},
                                                   seqan3::align_cfg::free_end_gaps_sequence1_trailing{true},
                                                   seqan3::align_cfg::free_end_gaps_sequence2_trailing{false}}
                | seqan3::align_cfg::edit_scheme
                | seqan3::align_cfg::parallel{static_cast<uint32_t>(seqan3::contrib::bgzf_thread_count)};

    // build pairs of reference with each haplotype
    auto && seq = seqan3::views::zip(seqan3::views::repeat(reference), seqan3::views::elements<1>(haplotypes));
    auto hap_iter = haplotypes.begin();
    for (auto ali : seqan3::align_pairwise(seq, config))
    {
        if (ali.score() != 0)
        {
            int pos = region_start + ali.sequence1_begin_position();
            seqan3::dna5_vector bufR{};
            seqan3::dna5_vector bufH{};
            for (auto && [ref, hap] : seqan3::views::zip(std::get<0>(ali.alignment()), std::get<1>(ali.alignment())))
            {
                if (ref != hap)
                {
                    if (ref != seqan3::gap{})
                    {
                        bufR.push_back(seqan3::dna5{}.assign_char(ref.to_char()));
                        ++pos;
                    }
                    if (hap != seqan3::gap{})
                        bufH.push_back(seqan3::dna5{}.assign_char(hap.to_char()));
                }
                else
                {
                    store_snp(junctions, bufR, bufH, pos++, ref_name, hap_iter->first);
                }
            }
            store_snp(junctions, bufR, bufH, pos, ref_name, hap_iter->first);
        }
        ++hap_iter;
    }
}

void detect_snp_and_indel(std::set<Junction> & junctions,
                          std::map<std::string, int32_t> & references_lengths,
                          std::filesystem::path const & genome_filename,
                          std::filesystem::path const & reads_filename,
                          uint64_t min_var_length)
{
    // Read the genome sequences.
    Genome genome = read_genome(genome_filename);

    // Fill the dictionary with reference information.
    for (auto && [genseq, name] : seqan3::views::zip(genome.seqs, genome.names))
    {
        auto && [iter, succ] = references_lengths.try_emplace(name, static_cast<int>(genseq.size()));
        if (!succ && static_cast<size_t>(iter->second) != genseq.size())
            throw std::runtime_error("contradicting sequence lengths encountered");
    }

    // Create activity profile.
    Activity activity = analyze_activity(reads_filename, genome, min_var_length);

    // Extract active regions from activity profile.
    std::vector<ActiveRegions> regions = get_active_regions(activity);

    // Analyze each active region and calculate haplotype sequences.
    for (auto && [ref_id, regv, genseq, name] : seqan3::views::zip(activity.refmap, regions, genome.seqs, genome.names))
    {
        seqan3::debug_stream << "Active regions for " << name << ": " << regv << std::endl;
        for (auto const & region : regv)
        {
            // Cut the region from the genome and convert to dna4.
            auto && ref = genseq | seqan3::views::slice(region.first, region.second)
                                 | seqan3::views::convert<seqan3::dna4>
                                 | seqan3::ranges::to<std::vector>();

            // Obtain haplotype sequences.
            auto haplotypes = find_haplotypes(ref_id, ref, region, reads_filename);
            // Call SNPs from the haplotypes.
            call_snp(junctions, std::move(ref), std::move(haplotypes), region.first, name);
        }
    }
}
