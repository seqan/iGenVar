#include <seqan3/io/stream/concept.hpp>
#include <fstream>
#include <map>

class info_entry
{
public:
    info_entry(std::string keyi, std::uint8_t numberi, std::string typei, std::string descriptioni, std::string sourcei, std::string versioni) :
        key{std::move(keyi)}, number{numberi}, type{std::move(typei)}, description{std::move(descriptioni)}, source{std::move(sourcei)}, version{std::move(versioni)}
    {}
    std::string key{};
    std::uint8_t number{};
    std::string type{};
    std::string description{};
    std::string source{};
    std::string version{};
};

class variant_header
{
public:
    variant_header() {}

    /*! \brief Set the file format for this VCF file header.
     *
     *  \param fileformati The input fileformat.
     */
    void set_fileformat(std::string fileformati)
    {
        fileformat = fileformati;
    }

    /*! \brief Add header information for a given INFO field.
     *
     *  \param info_keyi The INFO key name.
     *  \param numberi The number of values this key can hold.
     *  \param typei The type of values this key holds.
     *  \param descriptioni The description of this INFO field..
     *  \param sourcei The source of the INFO field.
     *  \param versioni The version of the source.
     */
    void add_meta_info(std::string info_keyi, std::uint8_t numberi, std::string typei, std::string descriptioni, std::string sourcei, std::string versioni)
    {
        info.push_back(info_entry{info_keyi, numberi, typei, descriptioni, sourcei, versioni});
    }

    /*! \brief Prints the VCF header to a given output.
     *
     * \tparam stream_type A stream to print the output to.
     *
     * \param out_stream The output stream to print to.
     */
    template<typename stream_type>
    //!cond
        requires seqan3::output_stream<stream_type>
    //!endcond
    void print(stream_type & out_stream)
    {
        out_stream << "##fileformat=" << fileformat << '\n';
        out_stream << "##source=" << source << '\n';
        for (const auto& i : info)
        {
            out_stream << "##INFO=<ID=" << i.key << ",Number=" << std::to_string(i.number) << ",Type=" << i.type << ",Description=\"" << i.description << "\",Source=\"" << i.source << "\",Version=\"" << i.version << "\">" << '\n';
        }
        out_stream << "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << '\n';
    }

private:
    std::string fileformat{"VCFv4.3"};
    std::string source{"iGenVarCaller"};
    std::vector<info_entry> info{};
};

class variant_record
{
public:
    variant_record() {}

    /*! \brief Set the chromosome for a variant.
     *
     * \param chromi The chromosome value to use.
     */
    void set_chrom(std::string chromi)
    {
        chrom = chromi;
    }

    /*! \brief Set the pos for a variant.
     *
     * \param posi The pos value to use.
     */
    void set_pos(std::uint64_t posi)
    {
        pos = posi;
    }

    /*! \brief Set the id for a variant.
     *
     * \param idi The id value to use.
     */
    void set_id(std::string idi)
    {
        id = idi;
    }

    /*! \brief Set the ref for a variant.
     *
     * \param refi The ref value to use.
     */
    void set_ref(std::string refi)
    {
        ref = refi;
    }

    /*! \brief Set the alt for a variant.
     *
     * \param alti The alt value to use.
     */
    void set_alt(std::string alti)
    {
        alt = alti;
    }

    /*! \brief Set the qual for a variant.
     *
     * \param quali The qual value to use.
     */
    void set_qual(float quali)
    {
        qual = quali;
    }

    /*! \brief Set the filter for a variant.
     *
     * \param filteri The filter value to use.
     */
    void set_filter(std::string filteri)
    {
        filter = filteri;
    }

    /*! \brief Add an INFO entry for a variant.
     *
     * \param info_key The INFO key to use.
     * \param info_value The INFO value to set.
     */
    void add_info(std::string info_key, std::string info_value)
    {
        info.insert_or_assign(info_key, info_value);
    }

    /*! \brief Prints the variant to a given output in VCF format.
     *
     * \tparam stream_type A stream to print the output to.
     *
     * \param out_stream The output stream to print to.
     */
    template<typename stream_type>
    //!cond
        requires seqan3::output_stream<stream_type>
    //!endcond
    void print(stream_type & out_stream)
    {
        std::map<std::string, std::string>::iterator last{};
        out_stream << chrom << '\t';
        out_stream << pos << '\t';
        out_stream << id << '\t';
        out_stream << ref << '\t';
        out_stream << alt << '\t';
        out_stream << qual << '\t';
        out_stream << filter << '\t';
        for (std::map<std::string, std::string>::iterator it = info.begin(); it != info.end(); ++it)
        {
            if (std::next(it) == info.end())
            {
                last = it;
                break;
            }
            out_stream << (*it).first << "=" << (*it).second << ";";
        }
        out_stream << (*last).first << "=" << (*last).second << '\n';
    }

private:
    std::string chrom{};
    std::uint64_t pos{};
    std::string id{"."};
    std::string ref{"N"};
    std::string alt{};
    float qual{};
    std::string filter{"PASS"};
    std::map<std::string, std::string> info{};
};
