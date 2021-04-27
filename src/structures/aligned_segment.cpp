#include "structures/aligned_segment.hpp"

#include <seqan3/alphabet/cigar/cigar.hpp>

int32_t AlignedSegment::get_reference_start() const
{
    return pos;
}

int32_t AlignedSegment::get_reference_end() const
{
    int32_t current_pos = pos;
    for (auto [element_length, element_operation] : cig)
    {
        switch(element_operation.to_char())
        {
            case 'M':
            case 'D':
            case 'N':
            case 'X':
            case '=':
                current_pos += element_length;
                break;
            case 'I':
            case 'S':
            case 'H':
            case 'P': // do nothing
                break;
            default:
                std::cerr << "The default case was accidentally triggered by the following element, which is no part of"
                             " a CIGAR string: "
                          << element_operation.to_char() << '\n';
                break;
        }
    }
    return current_pos;
}

int32_t AlignedSegment::get_left_soft_clip() const
{
    int32_t left_soft_clip = 0;
    for (auto [element_length, element_operation] : cig)
    {
        if(element_operation.to_char() == 'S')
        {
            left_soft_clip += element_length;
        }
        else if(element_operation.to_char() == 'M' ||
                element_operation.to_char() == '=' ||
                element_operation.to_char() == 'X' ||
                element_operation.to_char() == 'I')
        {
            break;
        }
    }
    return left_soft_clip;
}

int32_t AlignedSegment::get_right_soft_clip() const
{
    int32_t right_soft_clip = 0;
    for (auto [element_length, element_operation] : std::views::reverse(cig))
    {
        if(element_operation.to_char() == 'S')
        {
            right_soft_clip += element_length;
        }
        else if(element_operation.to_char() == 'M' ||
                element_operation.to_char() == '=' ||
                element_operation.to_char() == 'X' ||
                element_operation.to_char() == 'I')
        {
            break;
        }
    }
    return right_soft_clip;
}

int32_t AlignedSegment::get_query_start() const
{
    if (orientation == strand::forward)
    {
        return get_left_soft_clip();
    }
    else
    {
        return get_right_soft_clip();
    }
}

int32_t AlignedSegment::get_query_length() const
{
    int32_t current_length = 0;
    for (auto [element_length, element_operation] : cig)
    {
        switch(element_operation.to_char())
        {
            case 'M':
            case 'S':
            case 'I':
            case 'X':
            case '=':
                current_length += element_length;
                break;
            case 'D':
            case 'N':
            case 'H':
            case 'P': // do nothing
                break;
            default:
                std::cerr << "The default case was accidentally triggered by the following element, which is no part of"
                             " a CIGAR string: "
                          << element_operation.to_char() << '\n';
                break;
        }
    }
    return current_length;
}

int32_t AlignedSegment::get_query_end() const
{
    if (orientation == strand::forward)
    {
        return get_query_length() - get_right_soft_clip();
    }
    else
    {
        return get_query_length() - get_left_soft_clip();
    }
}

bool operator<(AlignedSegment const & lhs, AlignedSegment const & rhs)
{
    return lhs.get_query_start() != rhs.get_query_start()
            ? lhs.get_query_start() < rhs.get_query_start()
            : lhs.get_query_end() != rhs.get_query_end()
                ? lhs.get_query_end() < rhs.get_query_end()
                : lhs.mapq < rhs.mapq;
}

bool operator==(AlignedSegment const & lhs, AlignedSegment const & rhs)
{
    return lhs.orientation == rhs.orientation &&
           lhs.ref_name == rhs.ref_name &&
           lhs.pos == rhs.pos &&
           lhs.mapq == rhs.mapq &&
           lhs.cig == rhs.cig;
}
