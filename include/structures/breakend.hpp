#pragma once

enum struct sequence_type : uint8_t
{
    reference,
    read
};

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
    sequence_type seq_type;

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
inline stream_t operator<<(stream_t && stream, Breakend const & b)
{
    stream << ((b.seq_type == sequence_type::reference) ? "Reference" : "Read") << '\t'
           << b.seq_name << '\t'
           << b.position  << '\t'
           << ((b.orientation == strand::forward) ? "Forward" : "Reverse");
    return stream;
}

inline bool operator<(const Breakend & lhs, const Breakend & rhs)
{
    return lhs.seq_type != rhs.seq_type
            ? lhs.seq_type < rhs.seq_type
            : lhs.seq_name != rhs.seq_name
                ? lhs.seq_name < rhs.seq_name
                : lhs.orientation != rhs.orientation
                    ? lhs.orientation < rhs.orientation
                    : lhs.position < rhs.position;
}

inline bool operator==(const Breakend & lhs, const Breakend & rhs)
{
    return (lhs.seq_name == rhs.seq_name) &&
           (lhs.position == rhs.position) &&
           (lhs.orientation == rhs.orientation) &&
           (lhs.seq_type == rhs.seq_type);
}
