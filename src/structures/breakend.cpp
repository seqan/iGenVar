#include "structures/breakend.hpp"

bool operator<(const Breakend & lhs, const Breakend & rhs)
{
    return lhs.seq_type != rhs.seq_type
            ? lhs.seq_type < rhs.seq_type
            : lhs.seq_name != rhs.seq_name
                ? lhs.seq_name < rhs.seq_name
                : lhs.orientation != rhs.orientation
                    ? lhs.orientation < rhs.orientation
                    : lhs.position < rhs.position;
}

bool operator==(const Breakend & lhs, const Breakend & rhs)
{
    return (lhs.seq_name == rhs.seq_name) &&
           (lhs.position == rhs.position) &&
           (lhs.orientation == rhs.orientation) &&
           (lhs.seq_type == rhs.seq_type);
}
