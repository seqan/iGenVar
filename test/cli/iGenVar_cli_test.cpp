#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "cli_test.hpp"
#include <fstream>
#include <sstream>

std::string const default_alignment_long_reads_file_path = "simulated.minimap2.hg19.coordsorted_cutoff.sam";
std::string const vcf_out_file_path = "variants_file_out.vcf";

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
};

std::string const help_page_part_2
{
    "\n"
    "VERSION\n"
    "    Last update: 30-03-2021\n"
    "    iGenVar version: 0.0.3\n"
    "    SeqAn version: 3.0.3\n"
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
    "    -m, --method (List of detection_methods)\n"
    "          Choose the detection method(s) to be used. Default:\n"
    "          [cigar_string,split_read,read_pairs,read_depth]. Value must be one\n"
    "          of (method name or number)\n"
    "          [0,cigar_string,1,split_read,2,read_pairs,3,read_depth].\n"
    "    -c, --clustering_method (clustering_methods)\n"
    "          Choose the clustering method to be used. Default:\n"
    "          hierarchical_clustering. Value must be one of (method name or\n"
    "          number)\n"
    "          [0,simple_clustering,1,hierarchical_clustering,2,self_balancing_binary_tree,3,candidate_selection_based_on_voting].\n"
    "    -r, --refinement_method (refinement_methods)\n"
    "          Choose the refinement method to be used. Default: no_refinement.\n"
    "          Value must be one of (method name or number)\n"
    "          [0,no_refinement,1,sViper_refinement_method,2,sVirl_refinement_method].\n"
    "    -l, --min_var_length (signed 32 bit integer)\n"
    "          Specify what should be the minimum length of your SVs to be\n"
    "          detected. This value needs to be non-negative. Default: 30.\n"
    "    -x, --max_var_length (signed 32 bit integer)\n"
    "          Specify what should be the maximum length of your SVs to be\n"
    "          detected. SVs larger than this threshold can still be output as\n"
    "          translocations. This value needs to be non-negative. Default:\n"
    "          1000000.\n"
    "    -t, --max_tol_inserted_length (signed 32 bit integer)\n"
    "          Specify what should be the longest tolerated inserted sequence at\n"
    "          sites of non-INS SVs. This value needs to be non-negative. Default:\n"
    "          5.\n"
    "    -p, --max_overlap (signed 32 bit integer)\n"
    "          Specify the maximum allowed overlap between two alignment segments.\n"
    "          This value needs to be non-negative. Default: 10.\n"
};

// std::string expected_res_default
// {
//     "chr21\t41972615\tForward\tchr21\t41972616\tForward\t1\t1681\n"
//     "chr22\t17458417\tForward\tchr21\t41972615\tForward\t1\t2\n"
//     "chr22\t17458418\tForward\tchr21\t41972616\tForward\t2\t0\n"
// }

std::string expected_res_default
{
    "##fileformat=VCFv4.3\n"
    "##source=iGenVarCaller\n"
    "##contig=<ID=chr21,length=46709983>\n"
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
    "chr21\t41972616\t.\tN\t<INS>\t1\tPASS\tEND=41972616;SVLEN=1681;SVTYPE=INS\n"
};

std::string expected_res_empty
{
    "##fileformat=VCFv4.3\n"
    "##source=iGenVarCaller\n"
    "##contig=<ID=chr1,length=368>\n"
    "##INFO=<ID=SVTYPE,Number=1,Type=String,Description=\"Type of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=SVLEN,Number=1,Type=Integer,Description=\"Length of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "##INFO=<ID=END,Number=1,Type=Integer,Description=\"End position of SV called.\",Source=\"iGenVarCaller\",Version=\"1.0\">\n"
    "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n"
};

std::string expected_err_default_no_err
{
    "Detect junctions in long reads...\n"
    "INS: chr21\t41972615\tForward\tchr21\t41972616\tForward\t1681\tm2257/8161/CCS\n"
    "BND: chr21\t41972615\tReverse\tchr22\t17458415\tReverse\t2\tm41327/11677/CCS\n"
    "BND: chr21\t41972616\tReverse\tchr22\t17458416\tReverse\t0\tm21263/13017/CCS\n"
    "BND: chr21\t41972616\tReverse\tchr22\t17458416\tReverse\t0\tm38637/7161/CCS\n"
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
    cli_test_result result = execute_app("iGenVar", "-v");
    std::string expected_err
    {
        "[Error] Unknown option -v. In case this is meant to be a non-option/argument/parameter, please specify "
        "the start of non-options with '--'. See -h/--help for program information.\n"
    };
    EXPECT_EQ(result.exit_code, 65280);
    EXPECT_EQ(result.out, std::string{});
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
    cli_test_result result = execute_app("iGenVar", "-m 0");
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
                                         "-l -30");
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
                                         "-x -30");
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
                                         "-t -30");
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
                                         "-p -30");
    std::string expected_err
    {
        "[Error] You gave a negative max_overlap parameter.\n"
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
        "INS: chr21\t41972615\tForward\tchr21\t41972616\tForward\t1681\tm2257/8161/CCS\n"
        "The read depth method for long reads is not yet implemented.\n"
        "BND: chr21\t41972615\tReverse\tchr22\t17458415\tReverse\t2\tm41327/11677/CCS\n"
        "The read depth method for long reads is not yet implemented.\n"
        "BND: chr21\t41972616\tReverse\tchr22\t17458416\tReverse\t0\tm21263/13017/CCS\n"
        "The read depth method for long reads is not yet implemented.\n"
        "BND: chr21\t41972616\tReverse\tchr22\t17458416\tReverse\t0\tm38637/7161/CCS\n"
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

TEST_F(iGenVar_cli_test, with_detection_method_arguments)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "-m 0 -m 1");
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected_res_default);
    EXPECT_EQ(result.err, expected_err_default_no_err);
}

TEST_F(iGenVar_cli_test, with_detection_method_duplicate_arguments)
{
    cli_test_result result = execute_app("iGenVar",
                                         "-j", data(default_alignment_long_reads_file_path),
                                         "-m 0 -m 0");
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
                                         "-m 0 -m 1 -c 0 -r 0");
    std::string expected_err
    {
        "Detect junctions in long reads...\n"
        "INS: chr21\t41972615\tForward\tchr21\t41972616\tForward\t1681\tm2257/8161/CCS\n"
        "BND: chr21\t41972615\tReverse\tchr22\t17458415\tReverse\t2\tm41327/11677/CCS\n"
        "BND: chr21\t41972616\tReverse\tchr22\t17458416\tReverse\t0\tm21263/13017/CCS\n"
        "BND: chr21\t41972616\tReverse\tchr22\t17458416\tReverse\t0\tm38637/7161/CCS\n"
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
                                         "-m 9");
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
                                         "-m 2");

    // Check the output of junctions:
    seqan3::debug_stream << "Check the output of junctions... " << '\n';
    EXPECT_EQ(result.out, expected_res_empty);
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
                                         "-l 8 -m 0 -m 1");

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
