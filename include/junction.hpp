
#include <utility>
#include "breakend.hpp"

class junction
{
private:
    breakend mate1{};
    breakend mate2{};
    std::string read_name{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr junction()                    = default; //!< Defaulted.
    junction(junction const &)              = default; //!< Defaulted.
    junction(junction &&)                   = default; //!< Defaulted.
    junction & operator=(junction const &)  = default; //!< Defaulted.
    junction & operator=(junction &&)       = default; //!< Defaulted.
    ~junction()                             = default; //!< Defaulted.
    //!\}

    junction(breakend mate1, breakend mate2, std::string read_name) : mate1{std::move(mate1)},
                                                                      mate2{std::move(mate2)},
                                                                      read_name{std::move(read_name)}
    {
        if ((mate2.seq_type < mate1.seq_type) ||
            (mate2.seq_type == mate1.seq_type && mate2.seq_name < mate1.seq_name) ||
            (mate2.seq_type == mate1.seq_type && mate2.seq_name == mate1.seq_name && mate2.position < mate1.position))
        {
            std::swap(this->mate1, this->mate2);
            this->mate1.flip_orientation();
            this->mate2.flip_orientation();
        }
    }
    breakend get_mate1() const
    {
        return mate1;
    }
    breakend get_mate2() const
    {
        return mate2;
    }

    std::string get_read_name() const
    {
        return read_name;
    }
};

template <typename stream_t>
inline stream_t operator<<(stream_t && stream, junction const & junc)
{
    stream << junc.get_mate1() << '\t'
           << junc.get_mate2() << '\t'
           << junc.get_read_name();
    return stream;
}

inline bool operator<(const junction & lhs, const junction & rhs)
{
    return lhs.get_mate1() < rhs.get_mate1()
            ? true
            : rhs.get_mate1() < lhs.get_mate1()
                ? false
                : lhs.get_mate2() < rhs.get_mate2();
}

/*! \brief A junction is equal to another, if their mates are equal to each other. The read_name and supporting_reads
 *         are allowed to be unequal, because more than one read could support the same junction.
 *
 * \param lhs   left side junction
 * \param rhs   right side junction
 */
inline bool operator==(const junction & lhs, const junction & rhs)
{
    return (lhs.get_mate1() == rhs.get_mate1()) && (lhs.get_mate2() == rhs.get_mate2());
}
