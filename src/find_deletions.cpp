#include <fstream>
#include <seqan3/argument_parser/all.hpp>   // includes all necessary headers
#include "junction.hpp"

using namespace std;

vector<junction> read_junctions(const filesystem::path & junction_file_path)
{
    fstream junction_file;
    junction_file.open(junction_file_path,ios::in);
    vector<junction> junctions;
    if (junction_file.is_open()){
        string line, word;
        while(getline(junction_file, line)){
            stringstream s(line);
            vector<string> fields;
            while (getline(s, word, '\t')) { 
                fields.push_back(word);
            }
            if (fields.size() == 9)
            {
                breakend mate1{fields[1],
                               stoi(fields[2]),
                               (fields[3] == "Forward") ? strand::forward : strand::reverse,
                               (fields[0] == "Reference") ? sequence_type::reference : sequence_type::read};
                breakend mate2{fields[5],
                               stoi(fields[6]),
                               (fields[7] == "Forward") ? strand::forward : strand::reverse,
                               (fields[4] == "Reference") ? sequence_type::reference : sequence_type::read};
                junction new_junction{move(mate1), move(mate2), fields[8]};
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

struct cmd_arguments
{
    filesystem::path junction_file_path{};
};

void initialize_argument_parser(seqan3::argument_parser & parser, cmd_arguments & args)
{
    parser.info.author = "David Heller";
    parser.info.short_description = "Find deletions";
    parser.info.version = "0.0.1";
    parser.add_positional_option(args.junction_file_path, "Input junctions tab-separated format.");
}

int main(int argc, char ** argv)
{
    seqan3::argument_parser myparser{"partitionJunctions", argc, argv};        // initialise myparser
    cmd_arguments args{};
    initialize_argument_parser(myparser, args);
    try
    {
        myparser.parse();                                          // trigger command line parsing
    }
    catch (seqan3::parser_invalid_argument const & ext)                     // catch user errors
    {
        seqan3::debug_stream << "[Error] " << ext.what() << "\n"; // customise your error message
        return -1;
    }
    vector<junction> junctions = read_junctions(args.junction_file_path);
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
    return 0;
}