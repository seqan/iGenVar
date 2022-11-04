#include <gtest/gtest.h>

#include <fstream>
#include <vector>

#include "variant_detection/snp_indel_detection.hpp"

TEST(activity_detection, input_empty)
{
    Activity activity{};
    auto regions = get_active_regions(activity);
    EXPECT_TRUE(regions.empty());
}

TEST(activity_detection, input_zero)
{
    Activity activity{};
    activity.values.emplace_back(150, 0U); // 150x 0
    auto regions = get_active_regions(activity);
    EXPECT_EQ(regions.size(), 1U);
    EXPECT_TRUE(regions.front().empty());
}

TEST(activity_detection, input_active_everywhere)
{
    Activity activity{};
    activity.values.emplace_back(150, 9999U); // 150x 9999
    auto regions = get_active_regions(activity);
    EXPECT_EQ(regions.size(), 1U);
    EXPECT_EQ(regions.front().size(), 1U);
    EXPECT_EQ(regions.front().front(), (std::pair<int, int>{0, 150}));
}

TEST(activity_detection, input_interval)
{
    Activity activity{};
    activity.values.emplace_back(150, 0U);
    for (size_t pos = 34U; pos <= 44U; ++pos)
        activity.values.front()[pos] = 6;
    for (size_t pos = 134U; pos <= 144U; ++pos)
        activity.values.front()[pos] = 1;
    auto regions = get_active_regions(activity, 5);
    EXPECT_EQ(regions.size(), 1U);
    EXPECT_EQ(regions.front().size(), 2U);
    EXPECT_EQ(regions.front()[0], (std::pair<int, int>{30, 50}));
    EXPECT_EQ(regions.front()[1], (std::pair<int, int>{131, 149}));
}

TEST(activity_analysis, incompatible_genome)
{
    Genome genome;
    std::filesystem::path const reads = DATADIR"simulated.minimap2.hg19.coordsorted_cutoff.sam";
    // genome has too few sequences (0 instead of 1)
    EXPECT_THROW(analyze_activity(reads, genome, 20UL), std::runtime_error);

    using namespace seqan3::literals;
    genome.names.emplace_back("chr8");
    genome.seqs.push_back("ACG"_dna5);
    // genome length mismatch (3 instead of 46709983)
    EXPECT_THROW(analyze_activity(reads, genome, 20UL), std::runtime_error);
}
#include <seqan3/core/debug_stream.hpp>
TEST(activity_analysis, skip_ref)
{
    // Create a SAM file with three references.
    std::filesystem::path reads{std::filesystem::temp_directory_path()/"reads.sam"};
    {
        std::ofstream samfile(reads);
        samfile << "@HD\tVN:1.6\tSO:coordinate\n"
                << "@SQ\tSN:chr1\tLN:4\n"
                << "@SQ\tSN:chr2\tLN:3\n"
                << "test1\t16\tchr2\t1\t60\t10M\t=\t1\t0\tG\tF\n"; // chr1 skipped
        samfile.close();
    }

    Genome genome;
    using namespace seqan3::literals;
    genome.names.emplace_back("chr1");
    genome.seqs.push_back("ACGT"_dna5);
    genome.names.emplace_back("chr2");
    genome.seqs.push_back("ACG"_dna5);

    auto act = analyze_activity(reads, genome, 20UL);
    EXPECT_EQ(act.refmap, (std::vector<int>{-1, 1})); // first unused, second maps to seq 1
    EXPECT_TRUE(act.values[0].empty()); // first is unused => empty
    EXPECT_EQ(act.values[1], (std::vector<unsigned>{0, 0, 0})); // second has length 3 (=ref len), no activity added
}

TEST(store_snp, same_junction)
{
    using namespace seqan3::literals;
    std::set<Junction> junctions;
    seqan3::dna5_vector bufR = "ACG"_dna5;
    seqan3::dna5_vector bufH = "CT"_dna5;
    int pos = 30;
    std::string const ref_name = "chr1";

    store_snp(junctions, bufR, bufH, pos, ref_name, 0.1);
    bufR = "ACG"_dna5;
    bufH = "CG"_dna5;
    store_snp(junctions, bufR, bufH, pos, ref_name, 0.6); // this junction compares smaller
    bufR = "ACG"_dna5;
    bufH = "CT"_dna5;
    store_snp(junctions, bufR, bufH, pos, ref_name, 0.3); // add the initial junction again

    EXPECT_EQ(junctions.size(), 2UL);
    auto jctIt = junctions.begin(); // check 1st junction
    EXPECT_FLOAT_EQ(jctIt->get_quality(), 0.6F);
    ++jctIt; // go to 2nd junction
    EXPECT_FLOAT_EQ(jctIt->get_quality(), 0.4F); // 0.1 + 0.3
    EXPECT_TRUE(*jctIt != *junctions.begin());
}
