#include "variant_detection/variant_output.hpp"

#include <iostream> // for std::cout

#include "structures/junction.hpp"              // for class Junction
#include "variant_parser/variant_record.hpp"    // for class variant_header

void find_and_output_variants(std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::ostream & out_stream)
{
    variant_header header{};
    header.set_fileformat("VCFv4.3");
    header.add_meta_info("SVTYPE", 1, "String", "Type of SV called.", "iGenVarCaller", "1.0");
    header.add_meta_info("SVLEN", 1, "Integer", "Length of SV called.", "iGenVarCaller", "1.0");
    header.add_meta_info("END", 1, "Integer", "End position of SV called.", "iGenVarCaller", "1.0");
    header.print(out_stream);
    for (size_t i = 0; i < clusters.size(); ++i)
    {
        Breakend mate1 = clusters[i].get_average_mate1();
        Breakend mate2 = clusters[i].get_average_mate2();
        if (mate1.orientation == mate2.orientation)
        {
            if (mate1.seq_name == mate2.seq_name)
            {
                int32_t mate1_pos = mate1.position;
                int32_t mate2_pos = mate2.position;
                int32_t insert_size = clusters[i].get_average_inserted_sequence_size();
                if (mate1.orientation == strand::forward)
                {
                    int32_t distance = mate2_pos - mate1_pos;
                    //Deletion
                    if (distance >= args.min_var_length &&
                        distance < 1000000 &&
                        insert_size < 5)
                    {
                        variant_record tmp{};
                        tmp.set_chrom(mate1.seq_name);
                        tmp.set_qual(60);
                        tmp.set_alt("<DEL>");
                        tmp.add_info("SVTYPE", "DEL");
                        tmp.set_pos(mate1_pos);
                        tmp.add_info("SVLEN", std::to_string(-distance));
                        tmp.add_info("END", std::to_string(mate2_pos));
                        tmp.print(out_stream);
                    }
                    //Insertion
                    else if (distance == 1 &&
                             insert_size >= args.min_var_length)
                    {
                        variant_record tmp{};
                        tmp.set_chrom(mate1.seq_name);
                        tmp.set_qual(60);
                        tmp.set_alt("<INS>");
                        tmp.add_info("SVTYPE", "INS");
                        tmp.set_pos(mate1_pos);
                        tmp.add_info("SVLEN", std::to_string(insert_size));
                        tmp.add_info("END", std::to_string(mate1_pos));
                        tmp.print(out_stream);
                    }
                }
            }
        }
    }
}

//!\overload
void find_and_output_variants(std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::filesystem::path const & output_file_path)
{
    if (output_file_path.empty())
    {
        find_and_output_variants(clusters, args, std::cout);
    }
    else
    {
        std::ofstream out_file{output_file_path.c_str()};
        if (!out_file.good() || !out_file.is_open())
        {
            throw std::runtime_error{"Could not open file '" + output_file_path.string() + "' for reading."};
        }
        find_and_output_variants(clusters, args, out_file);
        out_file.close();
    }
}
