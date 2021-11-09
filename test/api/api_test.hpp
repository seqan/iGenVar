#pragma once

#include <gtest/gtest.h>

#include <memory>   // for std::unique_ptr
#include <vector>   // for std::vector

#include <seqan3/core/debug_stream.hpp> // include for debugging

#include "iGenVar.hpp"              // for global variable gVerbose
#include "structures/junction.hpp"  // for class Junction

#include "bamit/all.hpp" // for bamit index

/* From a discussion we decided to add a scope guard:
 * https://github.com/seqan/iGenVar/pull/169#pullrequestreview-774811822
 * There might be subtle problem with gVerbose = true. By default it is initialized to false. Depending on the unit test
 * order, it might be set to true by a previous unittest. So it is necessary to reset it to false with the help of a
 * scope guard function.
 * Idea taken from here:
 * https://stackoverflow.com/questions/28729545/abusing-c11-unique-ptr-to-execute-code-upon-leaving-the-scope)
 */

// Uses a unique_ptr to implement a scope_guard for guarding the gVerbose state
auto verbose_guard(bool verbose)
{
    // Lambda is executed when the scope of the unique_ptr ends
    auto onDestruction = [lastState = gVerbose](void*) {
        gVerbose = lastState;
    };
    gVerbose = verbose;
    return std::unique_ptr<void, decltype(onDestruction)>{(void*)1, onDestruction};
}

// This is a helper function for debugging, it prints junctions for comparing
void print_compare_junction_vectors(std::vector<Junction> const & junctions_expected_res,
                                    std::vector<Junction> const & junctions_res)
{
    for (size_t i = 0; i < junctions_res.size(); ++i)
    {
        seqan3::debug_stream << "------------------------------------------------------------------------------\n    "
                             << junctions_expected_res[i].get_mate1() << "\n == " << junctions_res[i].get_mate1() << "\n    "
                             << junctions_expected_res[i].get_mate2() << "\n == " << junctions_res[i].get_mate2() << "\n"
                             << junctions_expected_res[i].get_inserted_sequence() << " == " << junctions_res[i].get_inserted_sequence() << "\n"
                             << junctions_expected_res[i].get_tandem_dup_count() << " == " << junctions_res[i].get_tandem_dup_count() << "\n";
    }
}

// Helper function for comparing two bamit trees.
void compare_bamit_trees(std::unique_ptr<bamit::IntervalNode> const & t1, std::unique_ptr<bamit::IntervalNode> const & t2)
{
    if (!t1 && !t2) return;
    if (!t1 || !t2) EXPECT_EQ(1, 0);

    EXPECT_EQ(std::make_tuple(t1->get_start(), t1->get_end(), t1->get_file_position()),
              std::make_tuple(t2->get_start(), t2->get_end(), t2->get_file_position()));
    compare_bamit_trees(t1->get_left_node(), t2->get_left_node());
    compare_bamit_trees(t1->get_right_node(), t2->get_right_node());
}
