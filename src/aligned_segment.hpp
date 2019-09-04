#include <seqan3/alphabet/cigar/cigar.hpp>

struct aligned_segment
{
    int32_t ref_id;
    int32_t pos;
    strand orientation;
    std::vector<seqan3::cigar> cig;
    int32_t mapq;

    int32_t get_reference_start() const
    {
        return pos;
    }

    int32_t get_reference_end() const
    {
        int32_t current_pos = pos;
        for (auto [element_length, element_operation] : cig)
        {
            switch(element_operation.to_char())
            {
                case 'M':
                    current_pos += element_length;
                    break;
                case 'D':
                    current_pos += element_length;
                    break;
                case 'N':
                    current_pos += element_length;
                    break;
                case 'X':
                    current_pos += element_length;
                    break;
                case '=':
                    current_pos += element_length;
                    break;
            }
        }
        return current_pos;
    }

    int32_t get_left_soft_clip() const
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

    int32_t get_right_soft_clip() const
    {
        int32_t right_soft_clip = 0;
        for (auto [element_length, element_operation] : std::view::reverse(cig))
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

    int32_t get_query_start() const
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

    int32_t get_query_length() const
    {
        int32_t current_length = 0;
        for (auto [element_length, element_operation] : cig)
        {
            switch(element_operation.to_char())
            {
                case 'M':
                    current_length += element_length;
                    break;
                case 'S':
                    current_length += element_length;
                    break;
                case 'I':
                    current_length += element_length;
                    break;
                case 'X':
                    current_length += element_length;
                    break;
                case '=':
                    current_length += element_length;
                    break;
            }
        }
        return current_length;
    }

    int32_t get_query_end() const
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
};

template <typename stream_t>
inline stream_t operator<<(stream_t && stream, aligned_segment const & a)
{
    stream << a.ref_id << ";"
           << a.get_reference_start() << "-" << a.get_reference_end() << ";"
           << a.get_query_start() << "-" << a.get_query_end() << ";"
           << ((a.orientation == strand::forward) ? "+" : "-") << ";"
           << a.mapq;
    return stream;
}

inline bool operator<(const aligned_segment & lhs, const aligned_segment & rhs)
{
    return lhs.get_query_start() != rhs.get_query_start()
            ? lhs.get_query_start() < rhs.get_query_start()
            : lhs.get_query_end() != rhs.get_query_end()
                ? lhs.get_query_end() < rhs.get_query_end()
                : lhs.mapq < rhs.mapq;
}