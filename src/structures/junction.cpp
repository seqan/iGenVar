#include "structures/junction.hpp"

Breakend Junction::get_mate1() const
{
    return mate1;
}

Breakend Junction::get_mate2() const
{
    return mate2;
}

seqan3::dna5_vector Junction::get_inserted_sequence() const
{
    return inserted_sequence;
}

size_t Junction::get_tandem_dup_count() const
{
    return tandem_dup_count;
}

std::string Junction::get_read_name() const
{
    return read_name;
}

bool operator<(Junction const & lhs, Junction const & rhs)
{
    return lhs.get_mate1() != rhs.get_mate1()
           ? lhs.get_mate1() < rhs.get_mate1()
           : lhs.get_mate2() != rhs.get_mate2()
             ? lhs.get_mate2() < rhs.get_mate2()
             : lhs.get_tandem_dup_count() != rhs.get_tandem_dup_count()
               ? lhs.get_tandem_dup_count() < rhs.get_tandem_dup_count()
               : lhs.get_inserted_sequence() < rhs.get_inserted_sequence();
}

bool operator==(Junction const & lhs, Junction const & rhs)
{
    return (lhs.get_mate1() == rhs.get_mate1()) &&
           (lhs.get_mate2() == rhs.get_mate2()) &&
           (lhs.get_inserted_sequence() == rhs.get_inserted_sequence() &&
           (lhs.get_tandem_dup_count() == rhs.get_tandem_dup_count()));
}

bool operator!=(Junction const & lhs, Junction const & rhs)
{
    return !(lhs == rhs);
}
