#include "variant_detection/variant_output.hpp"

#include <iostream> // for std::cout

#include "structures/junction.hpp"              // for class Junction
#include "variant_parser/variant_record.hpp"    // for class variant_header

void find_and_output_variants(std::map<std::string, int32_t> & references_lengths,
                              std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::ostream & out_stream)
{
    variant_header header{};
    header.set_fileformat("VCFv4.3");
    header.add_meta_info("SVTYPE", 1, "String", "Type of SV called.", "iGenVarCaller", "1.0");
    header.add_meta_info("SVLEN", 1, "Integer", "Length of SV called.", "iGenVarCaller", "1.0");
    header.add_meta_info("END", 1, "Integer", "End position of SV called.", "iGenVarCaller", "1.0");
    header.print(references_lengths, args.vcf_sample_name, out_stream);
    for (size_t i = 0; i < clusters.size(); ++i)
    {
        size_t cluster_size = clusters[i].get_cluster_size();
        if (cluster_size >= args.min_qual)
        {
            Breakend mate1 = clusters[i].get_average_mate1();
            Breakend mate2 = clusters[i].get_average_mate2();
            if (mate1.orientation == mate2.orientation)
            {
                if (mate1.seq_name == mate2.seq_name)
                {
                    size_t insert_size = clusters[i].get_average_inserted_sequence_size();
                    if (mate1.orientation == strand::forward)
                    {
                        int distance = mate2.position - mate1.position - 1;
                        int sv_length = insert_size - distance;
                        // The SVLEN is neither too short nor too long than specified by the user.
                        if (std::abs(sv_length) >= args.min_var_length &&
                            std::abs(sv_length) <= args.max_var_length)
                        {
                            // Tandem Duplication
                            if (clusters[i].get_common_tandem_dup_count() > 0)
                            {
                                variant_record tmp{};
                                tmp.set_chrom(mate1.seq_name);
                                tmp.set_qual(cluster_size);
                                tmp.set_alt("<DUP:TANDEM>");
                                tmp.add_info("SVTYPE", "DUP");
                                // Increment position by 1 because VCF is 1-based
                                tmp.set_pos(mate1.position + 1);
                                tmp.add_info("SVLEN", std::to_string(distance));
                                // Increment end by 1 because VCF is 1-based
                                tmp.add_info("END", std::to_string(mate2.position + 1));
                                tmp.print(out_stream);
                            }
                            // Deletion (sv_length is negative)
                            else if (sv_length < 0 &&
                                     insert_size <= args.max_tol_inserted_length)
                            {
                                variant_record tmp{};
                                tmp.set_chrom(mate1.seq_name);
                                tmp.set_qual(cluster_size);
                                tmp.set_alt("<DEL>");
                                tmp.add_info("SVTYPE", "DEL");
                                // Increment position by 1 because VCF is 1-based
                                tmp.set_pos(mate1.position + 1);
                                tmp.add_info("SVLEN", std::to_string(sv_length));
                                // Increment end by 1 because VCF is 1-based
                                // Decrement end by 1 because deletion ends one base before mate2 begins
                                tmp.add_info("END", std::to_string(mate2.position));
                                tmp.print(out_stream);
                            }
                            // Insertion (sv_length is positive)
                            // for a small deletion inside of an insertion, the distance is a small negative value
                            else if (sv_length > 0 &&
                                     distance <= 0 &&
                                     std::abs(distance) <= args.max_tol_deleted_length)
                            {
                                variant_record tmp{};
                                tmp.set_chrom(mate1.seq_name);
                                tmp.set_qual(cluster_size);
                                tmp.set_alt("<INS>");
                                tmp.add_info("SVTYPE", "INS");
                                // Increment position by 1 because VCF is 1-based
                                tmp.set_pos(mate1.position + 1);
                                tmp.add_info("SVLEN", std::to_string(sv_length));
                                // Increment end by 1 because VCF is 1-based
                                tmp.add_info("END", std::to_string(mate1.position + 1));
                                tmp.print(out_stream);
                            }
                        }
                    }
                }
            }
        }
    }
}

//!\overload
void find_and_output_variants(std::map<std::string, int32_t> & references_lengths,
                              std::vector<Cluster> const & clusters,
                              cmd_arguments const & args,
                              std::filesystem::path const & output_file_path)
{
    if (output_file_path.empty())
    {
        find_and_output_variants(references_lengths, clusters, args, std::cout);
    }
    else
    {
        std::ofstream out_file{output_file_path.c_str()};
        if (!out_file.good() || !out_file.is_open())
        {
            throw std::runtime_error{"Could not open file '" + output_file_path.string() + "' for reading."};
        }
        find_and_output_variants(references_lengths, clusters, args, out_file);
        out_file.close();
    }
}
