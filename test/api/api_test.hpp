#pragma once

#include <gtest/gtest.h>

#include <memory>   // for std::unique_ptr

#include <seqan3/core/debug_stream.hpp> // include for debugging

#include "iGenVar.hpp"  // for global variable gVerbose

/* From a discussion we decided to add a scope guard:
 * https://github.com/seqan/iGenVar/pull/169#pullrequestreview-774811822
 * There might be subtle problem with gVerbose = true. By default it is initialized to false. Depending on the unit test
 * order, it might be set to true by a previous unittest. So it is necessary to reset it to false with the help of a
 * scope guard function.
 * Idea taken from here:
 * https://stackoverflow.com/questions/28729545/abusing-c11-unique-ptr-to-execute-code-upon-leaving-the-scope)
 */

// Uses a unique_ptr to implement a scope_guard for guarding the gVerbose state
auto verbose_guard(bool verbose) {
    // Lambda is being executed when the scope of the unique_ptr is being left
    auto onDestruction = [lastState = gVerbose](void*) {
        gVerbose = lastState;
    };
    gVerbose = verbose;
    return std::unique_ptr<void, decltype(onDestruction)>{(void*)1, onDestruction};
}
