#pragma once

#include <seqan3/argument_parser/auxiliary.hpp>     // for enumeration_names
#include <seqan3/argument_parser/validators.hpp>    // for value_list_validator

template <typename option_value_t>
class EnumValidator : public seqan3::value_list_validator<option_value_t>
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
     * \tparam range_type The type of range; must model std::ranges::forward_range and EnumValidator::option_value_type
     *                    must be constructible from the rvalue reference type of the given range.
     * \param[in] rng The range of valid values to test.
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

    //!\brief Returns a message that can be appended to the (positional) options help page info. <- todo...
    std::string get_help_page_message() const
    {
        std::vector<std::string> possible_values;
        for (auto && entry : values)
        {
            for (auto & [key, value] : seqan3::enumeration_names<option_value_t>)
            {
                if (entry == value)
                {
                    possible_values.emplace_back(key);
                }
            }
        }
        return seqan3::detail::to_string("Value must be one of (method name or number) ", possible_values, ".");
    }

private:

    //!\brief Minimum of the range to test.
    std::vector<option_value_type> values{};
};
