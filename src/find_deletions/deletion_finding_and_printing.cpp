#include "find_deletions/deletion_finding_and_printing.hpp"

#include <fstream>
/*! \brief Reads the input junction file and stores the junctions in a vector.
 *
 * \param junction_file_path input junction file
 *
 * \returns a vector of junctions
 */
std::vector<junction> read_junctions(std::filesystem::path const & junction_file_path)
{
    std::fstream junction_file;
    junction_file.open(junction_file_path, std::ios::in);
    std::vector<junction> junctions;
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
                breakend mate1{fields[1],
                               std::stoi(fields[2]),
                               (fields[3] == "Forward") ? strand::forward : strand::reverse,
                               (fields[0] == "Reference") ? sequence_type::reference : sequence_type::read};
                breakend mate2{fields[5],
                               std::stoi(fields[6]),
                               (fields[7] == "Forward") ? strand::forward : strand::reverse,
                               (fields[4] == "Reference") ? sequence_type::reference : sequence_type::read};
                junction new_junction{std::move(mate1), std::move(mate2), fields[8]};
                junctions.push_back(new_junction);
            }
        }
        junction_file.close();
    }
    return junctions;
}

/*! \brief Prints the header of a vcf file to a given outputfile.
 *
 * \param out_stream output stream object
 */
void print_vcf_header(std::ostream & out_stream)
{
    out_stream << "##fileformat=VCFv4.2" << '\n';
    out_stream << "##source=iGenVarCaller" << '\n';
    out_stream << "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << '\n';
}

/*! \brief Prints a deletion in vcf to a given outputfile.
 *
 * \param chrom chromosome, where the deletion is located
 * \param start start coordinate of the deletion
 * \param end   end coordinate of the deletion
 * \param qual  quality of the deletion, currently set to 60. ToDo: Requires a well-founded definition.
 * \param out_stream output stream object
 */
void print_deletion(std::string chrom, int32_t start, int32_t end, int32_t qual, std::ostream & out_stream)
{
    out_stream << chrom << '\t';
    out_stream << start << '\t';
    out_stream << "." << '\t';
    out_stream << "N" << '\t';
    out_stream << "<DEL>" << '\t';
    out_stream << qual << '\t';
    out_stream << "PASS" << '\t';
    out_stream << "SVTYPE=DEL;SVLEN=-" << end - start << ";END=" << end << '\n';
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
     std::vector<junction> junctions = read_junctions(junction_file_path);

     print_vcf_header(out_stream);
     for (size_t i = 0; i<junctions.size(); i++)
     {
         if (junctions[i].get_mate1().seq_type == sequence_type::reference &&
             junctions[i].get_mate2().seq_type == sequence_type::reference)
         {
             if (junctions[i].get_mate1().orientation == junctions[i].get_mate2().orientation)
             {
                 if (junctions[i].get_mate1().seq_name == junctions[i].get_mate2().seq_name)
                 {
                     int32_t mate1_pos = junctions[i].get_mate1().position;
                     int32_t mate2_pos = junctions[i].get_mate2().position;
                     if (junctions[i].get_mate1().orientation == strand::forward)
                     {
                         int32_t distance = mate2_pos - mate1_pos;
                         if (distance > 40 && distance < 100000)
                         {
                             print_deletion(junctions[i].get_mate1().seq_name, mate1_pos, mate2_pos, 60, out_stream);
                         }
                     }
                     else if (junctions[i].get_mate1().orientation == strand::reverse)
                     {
                         int32_t distance = mate1_pos - mate2_pos;
                         if (distance > 40 && distance < 100000)
                         {
                             print_deletion(junctions[i].get_mate1().seq_name, mate2_pos, mate1_pos, 60, out_stream);
                         }
                     }
                 }
             }
         }
     }
 }

//!\overload
void find_and_print_deletions(std::filesystem::path const & junction_file_path, std::filesystem::path const & output_file_path)
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
