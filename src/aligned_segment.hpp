#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

class AlignedSegment {
    int32_t ref_id;
    int32_t pos;
    bool strand;
    std::vector<cigar> cig;
    int32_t mapq;
    public:
        AlignedSegment (int32_t ref_id, int32_t pos, bool strand, std::vector<cigar> cig, int32_t mapq);
};

AlignedSegment::AlignedSegment (int32_t r, int32_t p, bool s, std::vector<cigar> c, int32_t m)
{
    ref_id = r;
    pos = p;
    strand = s;
    cig = c;
    mapq = m;
}