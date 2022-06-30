#include <gtest/gtest.h>

#include <fstream>
#include <iterator>

#include "structures/debruijn_graph.hpp"

TEST(debruijn_graph, empty)
{
    DeBruijnGraph graph;
    EXPECT_FALSE(graph.is_viable());
}

TEST(debruijn_graph, init)
{
    using namespace seqan3::literals;
    DeBruijnGraph graph;

    // the k-mer CG occurs twice
    EXPECT_FALSE(graph.init_sequence(2, "CGTCG"_dna4));

    // sequence as long as k
    EXPECT_TRUE(graph.init_sequence(5, "CGTCG"_dna4));
    EXPECT_FALSE(graph.is_viable());

    // sequence longer than k
    EXPECT_TRUE(graph.init_sequence(10, "CGTCGAAAATCAT"_dna4));
    EXPECT_FALSE(graph.is_viable());
}

TEST(debruijn_graph, haplotypes)
{
    using namespace seqan3::literals;
    DeBruijnGraph graph;
    EXPECT_TRUE(graph.init_sequence(4, "AAAACCCCGGGGUUUU"_dna4));
    EXPECT_FALSE(graph.is_viable());

    // add reads
    graph.add_read("AAAACCCCUUUU"_dna4);
    graph.add_read("AAAACCCCUUUU"_dna4);
    graph.add_read("AAAAGGGGUUUU"_dna4);
    graph.add_read("AAAAGGGGUUUU"_dna4);
    // not enough unique reads
    EXPECT_FALSE(graph.is_viable());
    graph.add_read("AUGC"_dna4); // add another unique k-mer
    graph.add_read("AUCG"_dna4); // add another unique k-mer
    EXPECT_TRUE(graph.is_viable());

    // 5 k-mers occur only once, prune them: CCCG CCGG CGGG AUGC AUCG
    graph.prune(2);
    auto haplotypes = graph.collect_haplotype_sequences();
    EXPECT_EQ(haplotypes.size(), 2U);

    EXPECT_FLOAT_EQ(haplotypes.front().first, .5F);
    EXPECT_EQ(haplotypes.front().second, "AAAAGGGGUUUU"_dna4);
    EXPECT_FLOAT_EQ(haplotypes.back().first, .5F);
    EXPECT_EQ(haplotypes.back().second, "AAAACCCCUUUU"_dna4);

    // export
    EXPECT_NO_THROW(graph.export_dot_format("final_graph.gv"));
    std::ifstream ifs("final_graph.gv");
    std::string file((std::istreambuf_iterator<char>(ifs)), std::istreambuf_iterator<char>());
    std::string expected = "digraph G {\n"
                           "42 [label=\"AGGG\"];\n"
                           "42 -> 170 [label=\"2\"];\n"
                           "10 [label=\"AAGG\"];\n"
                           "10 -> 42 [label=\"2\"];\n"
                           "2 [label=\"AAAG\"];\n"
                           "2 -> 10 [label=\"2\"];\n"
                           "127 [label=\"CTTT\"];\n"
                           "127 -> 255 [label=\"2\"];\n"
                           "95 [label=\"CCTT\"];\n"
                           "95 -> 127 [label=\"2\"];\n"
                           "87 [label=\"CCCT\"];\n"
                           "87 -> 95 [label=\"2\"];\n"
                           "255 [label=\"TTTT\"];\n"
                           "191 [label=\"GTTT\"];\n"
                           "191 -> 255 [label=\"2\"];\n"
                           "175 [label=\"GGTT\"];\n"
                           "175 -> 191 [label=\"2\"];\n"
                           "171 [label=\"GGGT\"];\n"
                           "171 -> 175 [label=\"2\"];\n"
                           "170 [label=\"GGGG\"];\n"
                           "170 -> 171 [label=\"2\"];\n"
                           "85 [label=\"CCCC\"];\n"
                           "85 -> 87 [label=\"2\"];\n"
                           "21 [label=\"ACCC\"];\n"
                           "21 -> 85 [label=\"2\"];\n"
                           "5 [label=\"AACC\"];\n"
                           "5 -> 21 [label=\"2\"];\n"
                           "1 [label=\"AAAC\"];\n"
                           "1 -> 5 [label=\"2\"];\n"
                           "0 [label=\"AAAA\"];\n"
                           "0 -> 2 [label=\"2\"];\n"
                           "0 -> 1 [label=\"2\"];\n"
                           "}\n";
    EXPECT_EQ(file, expected);
}
