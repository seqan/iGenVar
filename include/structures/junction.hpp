#pragma once

#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/alphabet/views/char_to.hpp>
#include <seqan3/alphabet/views/complement.hpp>
#include <seqan3/utility/views/to.hpp>

#include "structures/breakend.hpp"

class Junction
{
private:
    Breakend mate1{};
    Breakend mate2{};
    seqan3::dna5_vector inserted_sequence{};
    size_t tandem_dup_count{};
    std::string read_name{};

public:
    /*!\name Constructors, destructor and assignment
     * \{
     */
    constexpr Junction()                    = default; //!< Defaulted.
    Junction(Junction const &)              = default; //!< Defaulted.
    Junction(Junction &&)                   = default; //!< Defaulted.
    Junction & operator=(Junction const &)  = default; //!< Defaulted.
    Junction & operator=(Junction &&)       = default; //!< Defaulted.
    ~Junction()                             = default; //!< Defaulted.

    Junction(Breakend the_mate1,
             Breakend the_mate2,
             auto const & the_inserted_sequence,
             size_t the_tandem_dup_count,
             std::string the_read_name) : mate1{std::move(the_mate1)},
                                          mate2{std::move(the_mate2)},
                                          tandem_dup_count{the_tandem_dup_count},
                                          read_name{std::move(the_read_name)}
    {
        if ((mate2.seq_name < mate1.seq_name) ||
            (mate2.seq_name == mate1.seq_name && mate2.position < mate1.position))
        {
            std::swap(mate1, mate2);
            mate1.flip_orientation();
            mate2.flip_orientation();

            inserted_sequence = the_inserted_sequence | std::views::reverse | seqan3::views::complement | seqan3::views::to<seqan3::dna5_vector>;
        }
        else
        {
            inserted_sequence = the_inserted_sequence | seqan3::views::to<seqan3::dna5_vector>;
        }
    }
    //!\}

    //! \brief Returns the first mate of this junction.
    Breakend get_mate1() const;

    //! \brief Returns the second mate of this junction.
    Breakend get_mate2() const;

    /*! \brief Returns the sequence inserted between the two mates.
    *          If the two mates are connected directly, the inserted sequence is empty.
    */
    seqan3::dna5_vector get_inserted_sequence() const;

    //! \brief Returns the number of tandem copies of this junction.
    size_t get_tandem_dup_count() const;

    //! \brief Returns the name of the read giving rise to this junction.
    std::string get_read_name() const;
};

template <typename stream_t>
inline constexpr stream_t operator<<(stream_t && stream, Junction const & junc)
{
    stream << junc.get_mate1() << '\t'
           << junc.get_mate2() << '\t'
           << junc.get_inserted_sequence().size() << '\t'
           << junc.get_tandem_dup_count() << '\t'
           << junc.get_read_name();
    return stream;
}

/*! \brief A junction is smaller than another, if their first mate, second mate, or inserted sequence (in this order)
 *         is smaller than the corresponding element of the other junction.
 *
 * \param lhs   left side junction
 * \param rhs   right side junction
 */
bool operator<(Junction const & lhs, Junction const & rhs);

/*! \brief A junction is equal to another, if their mates and the inserted sequences are equal to each other.
 *         The read_name is allowed to be unequal, because more than one read could support the same junction.
 *
 * \param lhs - left side junction
 * \param rhs - right side junction
 */
bool operator==(Junction const & lhs, Junction const & rhs);

/*! \brief A junction is unequal to another, if their mates or the inserted sequences are unequal to each other.
 *
 * \param lhs - left side junction
 * \param rhs - right side junction
 */
bool operator!=(Junction const & lhs, Junction const & rhs);
