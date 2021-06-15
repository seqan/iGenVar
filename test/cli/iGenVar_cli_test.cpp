#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "cli_test.hpp"
#include <fstream>
#include <sstream>

std::string const default_alignment_long_reads_file_path = "simulated.minimap2.hg19.coordsorted_cutoff.sam";
std::string const vcf_out_file_path = "variants_file_out.vcf";
std::string const junctions_out_file_path = "junctions_file_out.txt";
std::string const clusters_out_file_path = "clusters_file_out.txt";

std::string const help_page_part_1
{
    "iGenVar - Detect genomic variants in a read alignment file\n"
    "==========================================================\n"
    "\n"
    "OPTIONS\n"
    "\n"
    "  Basic options:\n"
    "    -h, --help\n"
    "          Prints the help page.\n"
    "    -hh, --advanced-help\n"
    "          Prints the help page including advanced options.\n"
    "    --version\n"
    "          Prints the version information.\n"
    "    --copyright\n"
    "          Prints the copyright/license information.\n"
    "    --export-help (std::string)\n"
    "          Export the help page information. Value must be one of [html, man].\n"
    "    --version-check (bool)\n"
    "          Whether to check for the newest app version. Default: true.\n"
    "    -i, --input_short_reads (std::filesystem::path)\n"
    "          Input short read alignments in SAM or BAM format (Illumina).\n"
    "          Default: \"\". The input file must exist and read permissions must be\n"
    "          granted. Valid file extensions are: [sam, bam].\n"
    "    -j, --input_long_reads (std::filesystem::path)\n"
    "          Input long read alignments in SAM or BAM format (PacBio, Oxford\n"
    "          Nanopore, ...). Default: \"\". The input file must exist and read\n"
    "          permissions must be granted. Valid file extensions are: [sam, bam].\n"
    "    -o, --output (std::filesystem::path)\n"
    "          The path of the vcf output file. If no path is given, will output to\n"
    "          standard output. Default: \"\". Write permissions must be granted.\n"
    "          Valid file extensions are: [vcf].\n"
    "    -s, --vcf_sample_name (std::string)\n"
    "          Specify your sample name for the vcf header line. Default: MYSAMPLE.\n"
    "    -t, --threads (signed 16 bit integer)\n"
    "          Specify the number of decompression threads used for reading BAM\n"
    "          files. Default: 1.\n"
    "    -v, --verbose\n"
    "          If you set this flag, we provide additional details about what\n"
    "          iGenVar does. The detailed output is printed in the standard error.\n"
};

std::string const help_page_part_2
{
    "\n"
    "VERSION\n"
    "    Last update: 30-03-2021\n"
    "    iGenVar version: 0.0.3\n"
    "    SeqAn version: 3.1.0-rc.1\n"
    "\n"
    "URL\n"
    "    https://github.com/seqan/iGenVar/\n"
    "\n"
    "LEGAL\n"
    "    iGenVar Copyright: short_copyright\n"
    "    Author: Lydia Buntrock, David Heller, Joshua Kim\n"
    "    Contact: lydia.buntrock@fu-berlin.de\n"
    "    SeqAn Copyright: 2006-2021 Knut Reinert, FU-Berlin; released under the\n"
    "    3-clause BSDL.\n"
    "    For full copyright and/or warranty information see --copyright.\n"
};

std::string const help_page_advanced
{
    "    -a, --junctions (std::filesystem::path)\n"
    "          The path of the optional junction output file. If no path is given,\n"
    "          junctions will not be output. Default: \"\". Write permissions must be\n"
    "          granted.\n"
    "    -b, --clusters (std::filesystem::path)\n"
    "          The path of the optional cluster output file. If no path is given,\n"
    "          clusters will not be output. Default: \"\". Write permissions must be\n"
    "          granted.\n"
    "    -d, --method (List of detection_methods)\n"
    "          Choose the detection method(s) to be used. Value must be one of\n"
    "          (method name or number)\n"
    "          [0,cigar_string,1,split_read,2,read_pairs,3,read_depth]. Default:\n"
    "          [cigar_string,split_read,read_pairs,read_depth].\n"
    "    -c, --clustering_method (clustering_methods)\n"
    "          Choose the clustering method to be used. Value must be one of\n"
    "          (method name or number)\n"
    "          [0,simple_clustering,1,hierarchical_clustering,2,self_balancing_binary_tree,3,candidate_selection_based_on_voting].\n"
    "          Default: hierarchical_clustering.\n"
    "    -r, --refinement_method (refinement_methods)\n"
    "          Choose the refinement method to be used. Value must be one of\n"
    "          (method name or number)\n"
    "          [0,no_refinement,1,sViper_refinement_method,2,sVirl_refinement_method].\n"
    "          Default: no_refinement.\n"
    "    -k, --min_var_length (signed 32 bit integer)\n"
    "          Specify what should be the minimum length of your SVs to be\n"
    "          detected. This value needs to be non-negative. Default: 30.\n"
    "    -l, --max_var_length (signed 32 bit integer)\n"
    "          Specify what should be the maximum length of your SVs to be\n"
    "          detected. SVs larger than this threshold can still be output as\n"
    "          translocations. This value needs to be non-negative. Default:\n"
    "          100000.\n"
    "    -m, --max_tol_inserted_length (signed 32 bit integer)\n"
    "          Specify what should be the longest tolerated inserted sequence at\n"
    "          sites of non-INS SVs. This value needs to be non-negative. Default:\n"
    "          5.\n"
    "    -n, --max_overlap (signed 32 bit integer)\n"
    "          Specify the maximum allowed overlap between two alignment segments.\n"
    "          This value needs to be non-negative. Default: 10.\n"
    "    -q, --min_qual (signed 32 bit integer)\n"
    "          Specify the minimum quality (amount of supporting reads) of a\n"
    "          structural variant to be reported in the vcf output file. This value\n"
    "          needs to be non-negative. Default: 1.\n"
    "    -p, --partition_max_distance (signed 32 bit integer)\n"
    "          Specify the maximum distance in bp between members of the same\n"
    "          partition.This value needs to be non-negative. Default: 1000.\n"
    "    -w, --hierarchical_clustering_cutoff (double)\n"
    "          Specify the distance cutoff for the hierarchical clustering. This\n"
    "          value needs to be non-negative. Default: 0.5.\n"
};

// std::string expected_res_default
// {
//     "chr21\t41972615\tForward\tchr21\t41972616\tForward\t1\t1681\n"
//     "chr22\t17458417\tForward\tchr21\t41972615\tForward\t1\t2\n"
//     "chr22\t17458418\tForward\tchr21\t41972616\tForward\t2\t0\n"
// }

std::string expected_res_default_1
{
    "##fileformat=VCFv4.3\n"
    "##source=iGenVarCaller\n"
    "##contig=<ID=chr21,length=46709983>\n"
};

std::string general_header_lines
{
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##ALT=<ID=DEL,Number=1,Description=\"Deletion\">\n"
    "##ALT=<ID=DUP:TANDEM,Number=1,Description=\"Tandem Duplication\">\n"
    "##ALT=<ID=INS,Number=1,Description=\"Insertion of novel sequence\">\n"
    "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tMYSAMPLE\n"
};

std::string expected_res_default_2
{
    "chr21\t41972616\t.\tN\t<INS>\t1\tPASS\tEND=41972616;SVLEN=1681;SVTYPE=INS\tGT\t./.\n"
};

std::string expected_res_default = expected_res_default_1 + general_header_lines + expected_res_default_2;

std::string expected_res_empty
{
    "##fileformat=VCFv4.3\n"
    "##source=iGenVarCaller\n"
    "##contig=<ID=chr1,length=368>\n"
};

std::string expected_err_default_no_err
{
    "Detect junctions in long reads...\n"
    "Start clustering...\n"
    "Done with clustering. Found 2 junction clusters.\n"
    "No refinement was selected.\n"
};

TEST_F(iGenVar_cli_test, no_options)
{
    cli_test_result result = execute_app("iGenVar");
    std::string expected_res
    {
            "iGenVar - Detect genomic variants in a read alignment file\n"
            "==========================================================\n"
            "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected_res);
    EXPECT_EQ(result.err, std::string{});
}

// TODO (irallia): There is an open Issue, if we want to add the verbose option https://github.com/seqan/iGenVar/issues/20
TEST_F(iGenVar_cli_test, test_verbose_option)
{
    cli_test_result result = execute_app("iGenVar", "-j", data(default_alignment_long_reads_file_path), "--verbose");
    std::string expected_err
    {
        "Detect junctions in long reads...\n"
        "INS: chr21\t41972615\tForward\tchr21\t41972616\tForward\t1681\t0\tm2257/8161/CCS\n"
        "The read depth method for long reads is not yet implemented.\n"
        "BND: chr21\t41972615\tReverse\tchr22\t17458415\tReverse\t2\t0\tm41327/11677/CCS\n"
        "The read depth method for long reads is not yet implemented.\n"
        "BND: chr21\t41972616\tReverse\tchr22\t17458416\tReverse\t0\t0\tm21263/13017/CCS\n"
        "The read depth method for long reads is not yet implemented.\n"
        "BND: chr21\t41972616\tReverse\tchr22\t17458416\tReverse\t0\t0\tm38637/7161/CCS\n"
        "The read depth method for long reads is not yet implemented.\n"
        "Start clustering...\n"
        "Done with clustering. Found 2 junction clusters.\nNo refinement was selected.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected_res_default);
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, help_page_argument)
{
    cli_test_result result = execute_app("iGenVar", "-h");
    std::string expected_res = help_page_part_1 + help_page_part_2;

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected_res);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(iGenVar_cli_test, advanced_help_page_argument)
{
    cli_test_result result = execute_app("iGenVar", "-hh");
    std::string expected_res = help_page_part_1 + help_page_advanced + help_page_part_2;

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected_res);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(iGenVar_cli_test, fail_unknown_option)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "-y 0");

    std::string expected_err
    {
        "[Error] Unknown option -y. In case this is meant to be a non-option/argument/parameter, please specify the "
        "start of non-options with '--'. See -h/--help for program information.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, fail_missing_value)
{
    cli_test_result result = execute_app("iGenVar", "-i");
    std::string expected_err
    {
        "[Error] Missing value for option -i\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, fail_no_input_file)
{
    cli_test_result result = execute_app("iGenVar", "--method cigar_string");
    std::string expected_err
    {
        "[Error] You need to input at least one sam/bam file.\n"
        "Please use -i or -input_short_reads to pass a short read file or -j or -input_long_reads for a long read file.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, fail_negative_min_var_length)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--min_var_length -30");
    std::string expected_err
    {
        "[Error] You gave a negative min_var_length parameter.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, fail_negative_max_var_length)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--max_var_length -30");
    std::string expected_err
    {
        "[Error] You gave a negative max_var_length parameter.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, fail_negative_max_tol_inserted_length)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--max_tol_inserted_length -30");
    std::string expected_err
    {
        "[Error] You gave a negative max_tol_inserted_length parameter.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, fail_negative_max_overlap)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--max_overlap -30");
    std::string expected_err
    {
        "[Error] You gave a negative max_overlap parameter.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, fail_negative_min_qual)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--min_qual -30");
    std::string expected_err
    {
        "[Error] You gave a negative min_qual parameter.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, fail_negative_hierarchical_clustering_cutoff)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--hierarchical_clustering_cutoff -30");
    std::string expected_err
    {
        "[Error] You gave a negative hierarchical_clustering_cutoff parameter.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, with_default_arguments)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j ", data(default_alignment_long_reads_file_path));
    std::string expected_err
    {
        "Detect junctions in long reads...\n"
        "The read depth method for long reads is not yet implemented.\n"
        "The read depth method for long reads is not yet implemented.\n"
        "The read depth method for long reads is not yet implemented.\n"
        "The read depth method for long reads is not yet implemented.\n"
        "Start clustering...\n"
        "Done with clustering. Found 2 junction clusters.\n"
        "No refinement was selected.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected_res_default);
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, test_outfile)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j ", data(default_alignment_long_reads_file_path),
                                         "-o ", vcf_out_file_path);
    std::ifstream f;
    f.open(vcf_out_file_path);
    std::stringstream buffer;
    buffer << f.rdbuf();

    //this does not specifically check if file exists, rather if its readable.
    EXPECT_TRUE(f.is_open());
    EXPECT_EQ(buffer.str(), expected_res_default);
}

TEST_F(iGenVar_cli_test, test_intermediate_result_output)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j ", data(default_alignment_long_reads_file_path),
                                         "-a ", junctions_out_file_path,
                                         "-b ", clusters_out_file_path);
    std::ifstream f1;
    f1.open(junctions_out_file_path);
    std::stringstream buffer1;
    buffer1 << f1.rdbuf();

    // This does not specifically check if file exists, rather if its readable.
    EXPECT_TRUE(f1.is_open());
    EXPECT_NE(buffer1.str(), "");

    std::ifstream f2;
    f2.open(clusters_out_file_path);
    std::stringstream buffer2;
    buffer2 << f2.rdbuf();

    // This does not specifically check if file exists, rather if its readable.
    EXPECT_TRUE(f2.is_open());
    EXPECT_NE(buffer2.str(), "");
}

TEST_F(iGenVar_cli_test, with_detection_method_arguments)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--method cigar_string --method split_read");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected_res_default);
    EXPECT_EQ(result.err, expected_err_default_no_err);
}

TEST_F(iGenVar_cli_test, with_detection_method_duplicate_arguments)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--method cigar_string --method cigar_string");
    std::string expected_err
    {
        "[Error] The same detection method was selected multiple times.\n"
        "Methods to be used: [cigar_string,cigar_string]\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, test_direct_methods_input)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--method cigar_string --method split_read "
                                         "--clustering_method 0 --refinement_method 0");
    std::string expected_err
    {
        "Detect junctions in long reads...\n"
        "Start clustering...\n"
        "Done with clustering. Found 3 junction clusters.\n"
        "No refinement was selected.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected_res_default);
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, test_unknown_argument)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "--method 9");
    std::string expected_err
    {
        "[Error] You have chosen an invalid input value: 9. "
        "Please use one of: [0,cigar_string,1,split_read,2,read_pairs,3,read_depth]\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(iGenVar_cli_test, dataset_paired_end_mini_example)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-i", data("paired_end_mini_example.sam"),
                                         "--method read_pairs");

    // Check the output of junctions:
    seqan3::debug_stream << "Check the output of junctions... " << '\n';
    EXPECT_EQ(result.out, expected_res_empty + general_header_lines);
    seqan3::debug_stream << "done. " << '\n';

    // Check the debug output of junctions:
    seqan3::debug_stream << "Check the debug output of junctions... " << '\n';
    std::string expected_err
    {
        "Detect junctions in short reads...\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "The read pair method is not yet implemented.\n"
        "Start clustering...\n"
        "Done with clustering. Found 0 junction clusters.\n"
        "No refinement was selected.\n"
    };
    EXPECT_EQ(result.err, expected_err);
    seqan3::debug_stream << "done. " << '\n';
}

TEST_F(iGenVar_cli_test, dataset_single_end_mini_example)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data("single_end_mini_example.sam"),
                                         "--verbose",
                                         "--method cigar_string --method split_read",
                                         "--min_var_length 8 --max_var_length 400",
                                         "--hierarchical_clustering_cutoff 0.1");

    // Check the output of junctions:
    seqan3::debug_stream << "Check the output of junctions... " << '\n';
    std::ifstream output_res_file("../../data/output_res.txt");
    std::string output_res_str((std::istreambuf_iterator<char>(output_res_file)),
                                std::istreambuf_iterator<char>());
    EXPECT_EQ(result.out, output_res_str);
    seqan3::debug_stream << "done. " << '\n';

    // Check the debug output of junctions:
    seqan3::debug_stream << "Check the debug output of junctions... " << '\n';
    std::ifstream output_err_file("../../data/output_err.txt");
    std::string output_err_str((std::istreambuf_iterator<char>(output_err_file)),
                                std::istreambuf_iterator<char>());
    EXPECT_EQ(result.err, output_err_str);
    seqan3::debug_stream << "done. " << '\n';
}
