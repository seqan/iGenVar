#include "structures/junction.hpp"

Breakend Junction::get_mate1() const
{
    return mate1;
}

Breakend Junction::get_mate2() const
{
    return mate2;
}

std::string Junction::get_read_name() const
{
    return read_name;
}

bool operator<(const Junction & lhs, const Junction & rhs)
{
    return lhs.get_mate1() < rhs.get_mate1()
            ? true
            : rhs.get_mate1() < lhs.get_mate1()
                ? false
                : lhs.get_mate2() < rhs.get_mate2();
}

bool operator==(const Junction & lhs, const Junction & rhs)
{
    return (lhs.get_mate1() == rhs.get_mate1()) && (lhs.get_mate2() == rhs.get_mate2());
}
