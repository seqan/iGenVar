

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

struct breakend
{
    int32_t seq_id; // The id of the respective sequence
    int32_t position;
    strand orientation;
    sequence_type seq_type;
};

template <typename stream_t>
inline stream_t operator<<(stream_t && stream, breakend const & b)
{
    stream << "["
           << ((b.seq_type == sequence_type::reference) ? "Reference " : "Read ") << b.seq_id
           << "| Pos: " << b.position
           << "| Strand: " << ((b.orientation == strand::forward) ? "-->" : "<--")
           << "]";
    return stream;
}
