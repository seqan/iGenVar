#include <seqan3/alphabet/nucleotide/dna5.hpp>
#include <seqan3/io/sequence_file/input.hpp>
#include <seqan3/core/debug_stream.hpp>

#include "cli_test.hpp"
#include <fstream>
#include <sstream>

const std::string help_page_part_1
{
    "iGenVar - Detect junctions in a read alignment file\n"
    "===================================================\n"
    "\n"
    "POSITIONAL ARGUMENTS\n"
    "    ARGUMENT-1 (std::filesystem::path)\n"
    "          Input read alignments in SAM or BAM format. The input file must\n"
    "          exist and read permissions must be granted. Valid file extensions\n"
    "          are: [sam, bam].\n"
    "    ARGUMENT-2 (std::filesystem::path)\n"
    "          Output file for insertion alleles. Write permissions must be\n"
    "          granted. Valid file extensions are: [fa, fasta].\n"
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
};

const std::string help_page_part_2
{
    "\n"
    "VERSION\n"
    "    Last update: 19-01-2021\n"
    "    iGenVar version: 0.0.1\n"
    "    SeqAn version: 3.0.3\n"
    "\n"
    "URL\n"
    "    https://github.com/seqan/iGenVar/\n"
    "\n"
    "LEGAL\n"
    "    iGenVar Copyright: short_copyright\n"
    "    Author: David Heller & Lydia Buntrock\n"
    "    Contact: lydia.buntrock@fu-berlin.de\n"
    "    SeqAn Copyright: 2006-2021 Knut Reinert, FU-Berlin; released under the\n"
    "    3-clause BSDL.\n"
    "    For full copyright and/or warranty information see --copyright.\n"
};

const std::string help_page_advanced
{
    "    -m, --method (List of unsigned 8 bit integer)\n"
    "          Choose the detecting method(s) to be used. Default: [1,2,3,4]. Value\n"
    "          must be in range [1,4].\n"
    "    -c, --clustering_method (clustering_methods)\n"
    "          Choose the clustering method to be used. Default: simple_clustering.\n"
    "          Value must be one of (method name or number)\n"
    "          [simple_clustering,0,hierarchical_clustering,1,self_balancing_binary_tree,2,candidate_selection_based_on_voting,3].\n"
    "    -r, --refinement_method (refinement_methods)\n"
    "          Choose the refinement method to be used. Default: no_refinement.\n"
    "          Value must be one of (method name or number)\n"
    "          [no_refinement,0,sViper_refinement_method,1,2,sVirl_refinement_method].\n"
    "    -l, --min_var_length (unsigned 64 bit integer)\n"
    "          Specify what should be the minimum length of your SVs to be detected\n"
    "          (default 30 bp). Default: 30.\n"
};

TEST_F(detect_breakends, no_options)
{
    cli_test_result result = execute_app("detect_breakends");
    std::string expected
    {
            "iGenVar - Detect junctions in a read alignment file\n"
            "===================================================\n"
            "    Try -h or --help for more information.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(detect_breakends, fail_no_argument)
{
    cli_test_result result = execute_app("detect_breakends", "-v");
    std::string expected
    {
        "[Error] Unknown option -v. In case this is meant to be a non-option/argument/parameter, please specify "
        "the start of non-options with '--'. See -h/--help for program information.\n"
    };
    EXPECT_NE(result.exit_code, 0);
    EXPECT_EQ(result.out, std::string{});
    EXPECT_EQ(result.err, expected);
}

TEST_F(find_deletions, help_page_argument)
{
    cli_test_result result = execute_app("detect_breakends", "-h");
    std::string expected = help_page_part_1 + help_page_part_2;

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(find_deletions, advanced_help_page_argument)
{
    cli_test_result result = execute_app("detect_breakends", "-hh");
    std::string expected = help_page_part_1 + help_page_advanced + help_page_part_2;

    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, std::string{});
}

TEST_F(detect_breakends, with_arguments)
{
    cli_test_result result = execute_app("detect_breakends",
                                         data("simulated.minimap2.hg19.coordsorted_cutoff.sam"),
                                         "detect_breakends_insertion_file_out.fasta");
    std::string expected
    {
        "Reference\tchr21\t41972616\tForward\tRead\t0\t2294\tForward\t1\n"
        "Reference\tchr21\t41972616\tReverse\tRead\t0\t3975\tReverse\t1\n"
        "Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\t1\n"
        "Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\t2\n"
    };
    std::string expected_err
    {
        "INS1: Reference\tchr21\t41972616\tForward\tRead\t0\t2294\tForward\tm2257/8161/CCS\n"
        "INS2: Reference\tchr21\t41972616\tReverse\tRead\t0\t3975\tReverse\tm2257/8161/CCS\n"
        "The read pair method is not yet implemented.\n"
        "The read depth method is not yet implemented.\n"
        "BND: Reference\tchr22\t17458417\tForward\tReference\tchr21\t41972615\tForward\tm41327/11677/CCS\n"
        "The read pair method is not yet implemented.\n"
        "The read depth method is not yet implemented.\n"
        "BND: Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\tm21263/13017/CCS\n"
        "The read pair method is not yet implemented.\n"
        "The read depth method is not yet implemented.\n"
        "BND: Reference\tchr22\t17458418\tForward\tReference\tchr21\t41972616\tForward\tm38637/7161/CCS\n"
        "The read pair method is not yet implemented.\n"
        "The read depth method is not yet implemented.\n"
        "Start clustering...\n"
        "Done with clustering. Found 4 junction clusters.\n"
        "No refinement was selected.\n"
    };
    EXPECT_EQ(result.exit_code, 0);
    EXPECT_EQ(result.out, expected);
    EXPECT_EQ(result.err, expected_err);
}

TEST_F(detect_breakends, test_outfile)
{
    cli_test_result result = execute_app("detect_breakends",
                                         data("simulated.minimap2.hg19.coordsorted_cutoff.sam"),
                                         "detect_breakends_insertion_file_out.fasta");

    std::cout << "Current path is " << std::filesystem::current_path() << '\n';

    std::filesystem::path out_file_path = "detect_breakends_insertion_file_out.fasta";
    std::filesystem::path test_file_path = "../../data/detect_breakends_insertion_file_test.fasta";
    seqan3::sequence_file_input out_file{out_file_path};
    seqan3::sequence_file_input test_file{test_file_path};

    EXPECT_RANGE_EQ(out_file, test_file);
}
