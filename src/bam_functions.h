#ifndef INCLUDE_BAM_FUNCTIONS_H_
#define INCLUDE_BAM_FUNCTIONS_H_

enum BamFlags
{
    BAM_FLAG_MULTIPLE      = 0x0001,
    BAM_FLAG_ALL_PROPER    = 0x0002,
    BAM_FLAG_UNMAPPED      = 0x0004,
    BAM_FLAG_NEXT_UNMAPPED = 0x0008,
    BAM_FLAG_RC            = 0x0010,
    BAM_FLAG_NEXT_RC       = 0x0020,
    BAM_FLAG_FIRST         = 0x0040,
    BAM_FLAG_LAST          = 0x0080,
    BAM_FLAG_SECONDARY     = 0x0100,
    BAM_FLAG_QC_NO_PASS    = 0x0200,
    BAM_FLAG_DUPLICATE     = 0x0400,
    BAM_FLAG_SUPPLEMENTARY = 0x0800
};

inline bool
hasFlagUnmapped(uint32_t const & flag)
{
    return (flag & BAM_FLAG_UNMAPPED) == BAM_FLAG_UNMAPPED;
}

inline bool
hasFlagReverseComplement(uint32_t const & flag)
{
    return (flag & BAM_FLAG_RC) == BAM_FLAG_RC;
}

inline bool
hasFlagSecondary(uint32_t const & flag)
{
    return (flag & BAM_FLAG_SECONDARY) == BAM_FLAG_SECONDARY;
}

inline bool
hasFlagSupplementary(uint32_t const & flag)
{
    return (flag & BAM_FLAG_SUPPLEMENTARY) == BAM_FLAG_SUPPLEMENTARY;
}

inline bool
hasFlagDuplicate(uint32_t const & flag)
{
    return (flag & BAM_FLAG_DUPLICATE) == BAM_FLAG_DUPLICATE;
}

#endif  // #ifndef INCLUDE_BAM_FUNCTIONS_H_