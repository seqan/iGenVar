#include <gtest/gtest.h>

#include <vector>

#include "variant_detection/snp_indel_detection.hpp"

TEST(activity_analysis, input_empty)
{
    Activity activity{};
    auto regions = get_active_regions(activity);
    EXPECT_TRUE(regions.empty());
}

TEST(activity_analysis, input_zero)
{
    Activity activity{};
    activity.values.emplace_back(150, 0U); // 150x 0
    auto regions = get_active_regions(activity);
    EXPECT_EQ(regions.size(), 1U);
    EXPECT_TRUE(regions.front().empty());
}

TEST(activity_analysis, input_active_everywhere)
{
    Activity activity{};
    activity.values.emplace_back(150, 9999U); // 150x 9999
    auto regions = get_active_regions(activity);
    EXPECT_EQ(regions.size(), 1U);
    EXPECT_EQ(regions.front().size(), 1U);
    EXPECT_EQ(regions.front().front(), (std::pair<int, int>{0, 150}));
}

TEST(activity_analysis, input_interval)
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
