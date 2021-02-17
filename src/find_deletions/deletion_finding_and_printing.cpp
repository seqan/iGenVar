#include "find_deletions/deletion_finding_and_printing.hpp"

#include <iostream> // for std::cout

#include "structures/junction.hpp"              // for class Junction
#include "variant_parser/variant_record.hpp"    // for class variant_header

/*! \brief Reads the input junction file and stores the junctions in a vector.
 *
 * \param junction_file_path input junction file
 *
 * \returns a vector of junctions
 */
std::vector<Junction> read_junctions(std::filesystem::path const & junction_file_path)
{
    std::fstream junction_file;
    junction_file.open(junction_file_path, std::ios::in);
    std::vector<Junction> junctions;
    if (junction_file.is_open()){
        std::string line, word;
        while(getline(junction_file, line)){
            std::stringstream s(line);
            std::vector<std::string> fields;
            while (getline(s, word, '\t')) {
                fields.push_back(word);
            }
            if (fields.size() == 9)
            {
                Breakend mate1{fields[1],
                               std::stoi(fields[2]),
                               (fields[3] == "Forward") ? strand::forward : strand::reverse,
                               (fields[0] == "Reference") ? sequence_type::reference : sequence_type::read};
                Breakend mate2{fields[5],
                               std::stoi(fields[6]),
                               (fields[7] == "Forward") ? strand::forward : strand::reverse,
                               (fields[4] == "Reference") ? sequence_type::reference : sequence_type::read};
                Junction new_junction{std::move(mate1), std::move(mate2), fields[8]};
                junctions.push_back(new_junction);
            }
        }
        junction_file.close();
    }
    return junctions;
}

/*! \brief Detects deletions out of the junction file.

 * \cond
 * \param junction_file_path input junction file
 * \param out_stream output stream
 * \endcond
 *
 * \details Extracts deletions out of given breakends / junctions.
 */
void find_and_print_deletions(std::filesystem::path const & junction_file_path, std::ostream & out_stream)
{
    std::vector<Junction> junctions = read_junctions(junction_file_path);

    variant_header header{};
    header.set_fileformat("VCFv4.3");
    header.add_meta_info("SVTYPE", 1, "String", "Type of SV called.", "iGenVarCaller", "1.0");
    header.add_meta_info("SVLEN", 1, "Integer", "Length of SV called.", "iGenVarCaller", "1.0");
    header.add_meta_info("END", 1, "Integer", "End position of SV called.", "iGenVarCaller", "1.0");
    header.print(out_stream);
    for (size_t i = 0; i<junctions.size(); i++)
    {
        if (junctions[i].get_mate1().orientation == junctions[i].get_mate2().orientation)
        {
            if (junctions[i].get_mate1().seq_name == junctions[i].get_mate2().seq_name)
            {
                variant_record tmp{};
                tmp.set_chrom(junctions[i].get_mate1().seq_name);
                tmp.set_qual(60);
                tmp.set_alt("<DEL>");
                tmp.add_info("SVTYPE", "DEL");
                int32_t mate1_pos = junctions[i].get_mate1().position;
                int32_t mate2_pos = junctions[i].get_mate2().position;
                int32_t length{};
                if (junctions[i].get_mate1().orientation == strand::forward)
                {
                    int32_t distance = mate2_pos - mate1_pos;
                    if (distance > 40 && distance < 100000)
                    {
                        length = mate2_pos - mate1_pos;
                        tmp.set_pos(mate1_pos);
                        tmp.add_info("SVLEN", std::to_string(-length));
                        tmp.add_info("END", std::to_string(mate2_pos));
                        tmp.print(out_stream);
                        // print_deletion(junctions[i].get_mate1().seq_name, mate1_pos, mate2_pos, 60, out_stream);
                    }
                }
                else if (junctions[i].get_mate1().orientation == strand::reverse)
                {
                    int32_t distance = mate1_pos - mate2_pos;
                    if (distance > 40 && distance < 100000)
                    {
                        length = mate1_pos - mate2_pos;
                        tmp.set_pos(mate2_pos);
                        tmp.add_info("SVLEN", std::to_string(-length));
                        tmp.add_info("END", std::to_string(mate1_pos));
                        tmp.print(out_stream);
                        // print_deletion(junctions[i].get_mate1().seq_name, mate2_pos, mate1_pos, 60, out_stream);
                    }
                }
            }
        }
    }
}

//!\overload
void find_and_print_deletions(std::filesystem::path const & junction_file_path,
                              std::filesystem::path const & output_file_path)
{
    if (output_file_path.empty())
    {
        find_and_print_deletions(junction_file_path, std::cout);
    }
    else
    {
        std::ofstream out_file{output_file_path.c_str()};
        if (!out_file.good() || !out_file.is_open())
        {
            throw std::runtime_error{"Could not open file '" + output_file_path.string() + "' for reading."};
        }
        find_and_print_deletions(junction_file_path, out_file);
        out_file.close();
    }
}
