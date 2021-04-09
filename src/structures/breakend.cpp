#include "structures/breakend.hpp"

#include <tuple>

/*! \brief Compares two breakends.
 *
 * Breakends are compared in the following order:
 * 1. by the chromosome name
 * 2. by the orientation
 * 3. by their position
 */
bool operator<(Breakend const & lhs, Breakend const & rhs)
{
    return std::tie(lhs.seq_name, lhs.orientation, lhs.position) < std::tie(rhs.seq_name, rhs.orientation, rhs.position);
}

bool operator==(Breakend const & lhs, Breakend const & rhs)
{
    return (lhs.seq_name == rhs.seq_name) &&
           (lhs.position == rhs.position) &&
           (lhs.orientation == rhs.orientation);
}

bool operator!=(Breakend const & lhs, Breakend const & rhs)
{
    return !(lhs == rhs);
}
