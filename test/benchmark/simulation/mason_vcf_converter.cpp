#include <bio/alphabet/custom/char.hpp>
#include <bio/io/var/reader.hpp>
#include <bio/io/var/writer.hpp>

#include <seqan3/core/debug_stream.hpp> // for seqan3::debug_stream
#include <sharg/parser.hpp>             // for the SeqAn sharg::parser

#include <filesystem>                   // for filesystem
#include <variant>                      // for std::get_if

using namespace bio::alphabet::literals;
using namespace std::string_literals;

struct cmd_arguments
{
// Input:
    /* -i */ std::filesystem::path input_file_path{""};
// Output:
    /* -o */ std::filesystem::path output_file_path{};
};

/*! \brief Fill ALT fields in VCF with SV types.
 *
 * \param[in] args - command line arguments:\n
 *                   **args.input_file_path** - Mason VCF input file, path to the VCF file\n
 *                   **args.output_file_path** output file - path for the VCF file - *default: standard output*\n
 *
 * \details Takes SVTYPE of SV InDels and fills their ALT fields.
 */
void convert_vcf(cmd_arguments const & args)
{
    bio::io::var::reader_options opt{.record = bio::io::var::record_default{} }; // record_default ist deep; record_default_shallow ist shallow
    bio::io::var::reader reader{args.input_file_path, opt};
    bio::io::var::writer writer{args.output_file_path};

    bio::io::var::header hdr = reader.header();
    writer.set_header(hdr);
    hdr.formats.clear();
    hdr.formats.push_back({ .id = "GT", .number = 1, .type_id = bio::io::var::value_type_id::string,
                            .description = "Genotype"});

    for (auto & rec : reader) {
        rec.ref = "N"_dna5;
        rec.qual = 1;
        rec.genotypes.clear();
        rec.genotypes.push_back({.id = "GT", .value = std::vector{"./."s}});

        if (rec.id.starts_with("sim_sv_indel")) {
            // If we would assign a std::string to this within the inner if-statement, then that string goes out-of-scope
            // at the end of the if-statement, so the string_view is "dangling" when we push the record to the writer.
            std::string alt;
            int32_t end{-1};
            for (auto & info : rec.info) {
                if (info.id == "SVTYPE") {
                    if(std::string* s = std::get_if<std::string>(&info.value)) {
                        alt = "<" + *s + ">";
                        rec.alt = {alt};
                    }
                } else if (info.id == "SVLEN") {
                    if(auto* i = std::get_if<std::vector<int32_t>>(&info.value); i && i->size() == 1) {
                        if ((*i)[0] <= 0) // -> DEL
                            end = rec.pos - (*i)[0];
                        else // -> INS
                            end = rec.pos;

                    }
                }
            }
            rec.info.push_back({.id = "END", .value = end});
        }
        writer.push_back(rec);
    }
}

void initialize_argument_parser(sharg::parser & parser, cmd_arguments & args)
{
    parser.info.author = "Lydia Buntrock";
    parser.info.app_name = "MasonVcfConverter";
    parser.info.man_page_title = "Mason VCF converter";
    parser.info.short_description = "Fill ALT fields in VCF with SV types.";
    parser.info.version = "0.0.1";
    parser.info.date = "18-10-2022";    // last update
    parser.info.email = "lydia.buntrock@fu-berlin.de";
    parser.info.url = "https://github.com/seqan/iGenVar/";

    // Options - Input / Output:
    parser.add_option(args.input_file_path, sharg::config{
                        .short_id = 'i', .long_id = "input",
                        .description = "Input VCF file created by Mason.",
                        .validator = sharg::input_file_validator{{"vcf"}} });
    parser.add_option(args.output_file_path, sharg::config{
                        .short_id = 'o', .long_id = "output",
                        .description = "The path of the VCF output file. If no path is given, will output to standard output.",
                        .validator = sharg::output_file_validator{sharg::output_file_open_options::open_or_create, {"vcf"}}});
}

int main(int argc, char ** argv)
{
    sharg::parser parser{"MasonVcfConverter", argc, argv};          // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(parser, args);

    // Parse the given arguments and catch possible errors.
    try
    {
        parser.parse();                                             // trigger command line parsing
    }
    catch (sharg::parser_error const & ext)                         // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << '\n';   // customise your error message
        return -1;
    }

    // Check if we have an input file.
    if (args.input_file_path.empty())
    {
        seqan3::debug_stream << "[Error] You need to input a VCF file.\n"
                             << "Please use -i or -input to pass the VCF file.\n";
        return -1;
    }

    convert_vcf(args);

    return 0;
}
