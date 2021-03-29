#include "structures/breakend.hpp"

/*! \brief Compares two breakends.
 *
 * Breakends are compared in the following order:
 * 1. by the chromosome name
 * 2. by the orientation
 * 3. by their position
 */
bool operator<(Breakend const & lhs, Breakend const & rhs)
{
    return lhs.seq_name != rhs.seq_name
            ? lhs.seq_name < rhs.seq_name
            : lhs.orientation != rhs.orientation
                ? lhs.orientation < rhs.orientation
                : lhs.position < rhs.position;
}

bool operator==(Breakend const & lhs, Breakend const & rhs)
{
    return (lhs.seq_name == rhs.seq_name) &&
           (lhs.position == rhs.position) &&
           (lhs.orientation == rhs.orientation);
}
