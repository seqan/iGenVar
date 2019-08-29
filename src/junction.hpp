
#include "breakend.hpp"

class junction
{
public:

    junction(breakend mate1, breakend mate2) : mate1{std::move(mate1)}, mate2{std::move(mate2)}
    {}

private:
    breakend mate1{};
    breakend mate2{};
};


template <typename stream_t>
inline stream_t operator<<(stream_t && stream, junction const & junc)
{
    stream << junc.mate1 << " - " << junc.mate2;
    return stream;
}
