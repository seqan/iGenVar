#pragma once

#include <string>

enum struct strand : uint8_t
{
    forward,
    reverse
};

struct Breakend
{
    std::string seq_name; // The id of the respective sequence
    int32_t position;
    strand orientation;

    void flip_orientation()
    {
        if (orientation == strand::forward)
        {
            orientation = strand::reverse;
        }
        else
        {
            orientation = strand::forward;
        }
    }
};

template <typename stream_t>
inline constexpr stream_t operator<<(stream_t && stream, Breakend const & b)
{
    stream << b.seq_name << '\t'
           << b.position  << '\t'
           << ((b.orientation == strand::forward) ? "Forward" : "Reverse");
    return stream;
}

bool operator<(Breakend const & lhs, Breakend const & rhs);

bool operator==(Breakend const & lhs, Breakend const & rhs);

bool operator!=(Breakend const & lhs, Breakend const & rhs);
