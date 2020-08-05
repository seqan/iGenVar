#include "find_deletions/deletion_finding_and_printing.hpp"

std::vector<junction> read_junctions(const std::filesystem::path & junction_file_path)
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

void print_vcf_header()
{
    std::cout << "##fileformat=VCFv4.2" << '\n';
    std::cout << "##source=iGenVarCaller" << '\n';
    std::cout << "CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO" << '\n';
}

void print_deletion(std::string chrom, int32_t start, int32_t end, int32_t qual)
{
    std::cout << chrom << '\t';
    std::cout << start << '\t';
    std::cout << "." << '\t';
    std::cout << "N" << '\t';
    std::cout << "<DEL>" << '\t';
    std::cout << qual << '\t';
    std::cout << "PASS" << '\t';
    std::cout << "SVTYPE=DEL;SVLEN=-" << end - start << ";END=" << end << '\n';
}

void find_and_print_deletions(const std::filesystem::path & junction_file_path)
{
    std::vector<junction> junctions = read_junctions(junction_file_path);
    print_vcf_header();
    for (size_t i = 0; i<junctions.size(); i++)
    {
        if (junctions[i].get_mate1().seq_type == sequence_type::reference && junctions[i].get_mate2().seq_type == sequence_type::reference)
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
                            print_deletion(junctions[i].get_mate1().seq_name, mate1_pos, mate2_pos, 60);
                        }
                    }
                    else if (junctions[i].get_mate1().orientation == strand::reverse)
                    {
                        int32_t distance = mate1_pos - mate2_pos;
                        if (distance > 40 && distance < 100000)
                        {
                            print_deletion(junctions[i].get_mate1().seq_name, mate2_pos, mate1_pos, 60);
                        }
                    }
                }
            }
        }
    }
}
