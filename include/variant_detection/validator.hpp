#pragma once

#include <seqan3/std/algorithm>                             // for std::ranges::move

#include <seqan3/argument_parser/auxiliary.hpp>             // for enumeration_names
#include <seqan3/argument_parser/exceptions.hpp>            // for seqan3::validation_error
#include <seqan3/core/debug_stream/detail/to_string.hpp>    // for seqan3::detail::to_string

template <typename option_value_t>
class EnumValidator
{
public:
    //!\brief Type of values that are tested by validator
    using option_value_type = option_value_t;

    /*!\name Constructors, destructor and assignment
     * \{
     */
    EnumValidator() = default;                                  //!< Defaulted.
    EnumValidator(EnumValidator const &) = default;             //!< Defaulted.
    EnumValidator(EnumValidator &&) = default;                  //!< Defaulted.
    EnumValidator & operator=(EnumValidator const &) = default; //!< Defaulted.
    EnumValidator & operator=(EnumValidator &&) = default;      //!< Defaulted.
    ~EnumValidator() = default;                                 //!< Defaulted.

    /*!\brief Constructing from a range.
     * \tparam range_type - the type of range; must model std::ranges::forward_range and
     *                      EnumValidator::option_value_type must be constructible from the rvalue reference type of the
     *                      given range
     * \param[in] rng - the range of valid values to test
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires std::constructible_from<option_value_type, std::ranges::range_rvalue_reference_t<range_type>>
    //!\endcond
    EnumValidator(range_type rng)
    {
        values.clear();
        std::ranges::move(std::move(rng), std::cpp20::back_inserter(values));
        std::ranges::sort(values);
        values.erase(std::unique(values.begin(), values.end()), values.end());
    }
    //!\}

    /*!\brief Tests whether cmp lies inside values.
     * \param cmp The input value to check.
     * \throws seqan3::validation_error
     */
    void operator()(option_value_type const & cmp) const
    {
        if (!(std::find(values.begin(), values.end(), cmp) != values.end()))
            throw seqan3::validation_error{seqan3::detail::to_string("Value ", cmp, " is not one of ",
                                                                     std::views::all(values), ".")};
    }

    /*!\brief Tests whether every element in \p range lies inside values.
     * \tparam range_type The type of range to check; must model std::ranges::forward_range.
     * \param  range      The input range to iterate over and check every element.
     * \throws seqan3::validation_error
     */
    template <std::ranges::forward_range range_type>
    //!\cond
        requires std::convertible_to<std::ranges::range_value_t<range_type>, option_value_type>
    //!\endcond
    void operator()(range_type const & range) const
    {
        std::for_each(std::ranges::begin(range), std::ranges::end(range), [&] (auto cmp) { (*this)(cmp); });
    }

    //!\brief Returns a message that can be appended to the (positional) options help page info.
    std::string get_help_page_message() const
    {
        auto map = seqan3::enumeration_names<option_value_t>;
        std::vector<std::pair<std::string_view, option_value_t>> key_value_pairs(map.begin(), map.end());
        std::ranges::sort(key_value_pairs, [] (auto pair1, auto pair2)
            {
                if constexpr (std::totally_ordered<option_value_t>)
                {
                    if (pair1.second != pair2.second)
                        return pair1.second < pair2.second;
                }
                return pair1.first < pair2.first;
            });

        return seqan3::detail::to_string("Value must be one of (method name or number) ",
                                         key_value_pairs | std::views::keys,
                                         ".");
    }

private:

    //!\brief Minimum of the range to test.
    std::vector<option_value_type> values{};
};
