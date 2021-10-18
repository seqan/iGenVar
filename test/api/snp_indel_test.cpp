#include <gtest/gtest.h>

#include <vector>

#include "variant_detection/snp_indel_detection.hpp"

TEST(activity_analysis, input_empty)
{
    std::vector<unsigned> activity{};
    auto regions = active_regions(activity);
    EXPECT_TRUE(regions.empty());
}

TEST(activity_analysis, input_zero)
{
    std::vector<unsigned> activity(150, 0u);
    auto regions = active_regions(activity);
    EXPECT_TRUE(regions.empty());
}

TEST(activity_analysis, input_active_everywhere)
{
    std::vector<unsigned> activity(150, 9999u);
    auto regions = active_regions(activity);
    EXPECT_EQ(regions.size(), 1u);
    EXPECT_EQ(regions[0], (std::pair<size_t, size_t>{0, 149}));
}

TEST(activity_analysis, input_interval)
{
    std::vector<unsigned> activity(150, 0u);
    for (size_t pos = 34u; pos <= 44u; ++pos)
        activity[pos] = 6;
    for (size_t pos = 134u; pos <= 144u; ++pos)
        activity[pos] = 1;
    auto regions = active_regions(activity);
    EXPECT_EQ(regions.size(), 2u);
    EXPECT_EQ(regions[0], (std::pair<size_t, size_t>{30, 49}));
    EXPECT_EQ(regions[1], (std::pair<size_t, size_t>{131, 148}));
}
