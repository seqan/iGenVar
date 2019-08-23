#include <seqan3/core/debug_stream.hpp>

using namespace seqan3;

class Breakend {
    bool in_reference1;
    int32_t chromosome1;
    int32_t position1;
    bool forward1;
    bool in_reference2;
    int32_t chromosome2;
    int32_t position2;
    bool forward2;
    public:
        Breakend (bool, int32_t, int32_t, bool, bool, int32_t, int32_t, bool);
        void print_vcf_entry ();
} ;

Breakend::Breakend (bool in_ref1, int32_t chr1, int32_t pos1, bool fwd1, bool in_ref2, int32_t chr2, int32_t pos2, bool fwd2)
{
    in_reference1 = in_ref1;
    chromosome1 = chr1;
    position1 = pos1;
    forward1 = fwd1;
    in_reference2 = in_ref2;
    chromosome2 = chr2;
    position2 = pos2;
    forward2 = fwd2;
}

void Breakend::print_vcf_entry()
{
    printf("%d\t%d\t%s\t%d\t%d\t%s\n", chromosome1, position1, forward1 ? "-->" : "<--", chromosome2, position2, forward2 ? "-->" : "<--");
}