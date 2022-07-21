#include <algorithm>
#include <fstream>

#include <seqan3/search/views/kmer_hash.hpp>
#include <seqan3/utility/views/to.hpp>
#include <seqan3/utility/views/zip.hpp>

#include <lemon/connectivity.h>
#include <lemon/path.h>

#include "structures/debruijn_graph.hpp"

bool DeBruijnGraph::init_sequence(unsigned char kmer_len, seqan3::dna4_vector const & reference)
{
    // Initialize the graph.
    len = kmer_len;
    viable = false;
    graph.clear();
    assert(len <= reference.size());

    // Reserve memory for nodes and arcs, because we know the size a-priori.
    graph.reserveNode(static_cast<int>(reference.size() + 1 - len));
    graph.reserveArc(static_cast<int>(reference.size() - len));

    lemon::ListDigraph::Node prev_node = lemon::INVALID; // the previously processed node (initial state: none)
    for (size_t kmer : seqan3::views::kmer_hash(reference, seqan3::ungapped{len})) // iterate kmers
    {
        if (nodes(kmer) == lemon::INVALID) // the kmer is new
        {
            nodes.set(graph.addNode(), kmer); // add node
            if (prev_node != lemon::INVALID) // if there exists a previous node
                arcs.set(graph.addArc(prev_node, nodes(kmer)), 0); // add arc from previous to current node

            prev_node = nodes(kmer); // update the previously processed node
        }
        else // abort: the kmer already exists (we have found a cycle)
        {
            len = 0;
            return false;
        }
    }
    return true;
}

void DeBruijnGraph::add_read(seqan3::dna4_vector const & read)
{
    assert(len > 0); // adding reads to an uninitialized graph is an error
    lemon::ListDigraph::Node prev_node = lemon::INVALID; // the previously processed node
    for (size_t kmer : seqan3::views::kmer_hash(read, seqan3::ungapped{len})) // iterate kmers
    {
        if (nodes(kmer) == lemon::INVALID) // kmer is not in graph: new node
            nodes.set(graph.addNode(), kmer);

        if (prev_node != lemon::INVALID)
        {
            lemon::ListDigraph::Arc arc = lemon::findArc(graph, prev_node, nodes(kmer));
            if (arc != lemon::INVALID) // arc exists: increment read count
                ++arcs[arc];
            else
                arcs.set(graph.addArc(prev_node, nodes(kmer)), 1); // new arc
        }
        prev_node = nodes(kmer); // update processed node
    }
    viable = std::nullopt; // we do not know if graph is viable anymore
}

bool DeBruijnGraph::is_viable()
{
    if (!viable) // status is unknown
    {
        if (len == 0 || !lemon::dag(graph)) // graph is uninitialized or has cycles (following directed arcs)
        {
            viable = false;
        }
        else // check if graph is complex enough (number of non-unique nodes < 4 * number of unique nodes)
        {
            int num_unique{}; // number of unique kmers (0-1 in and out arcs)
            int num_non_unique{}; // number of non-unique kmers (>1 in or out arcs)
            for (lemon::ListDigraph::NodeIt node(graph); node != lemon::INVALID; ++node)
            {
                lemon::ListDigraph::InArcIt in_arc(graph, node);
                lemon::ListDigraph::OutArcIt out_arc(graph, node);
                // true if the node has either no ingoing arc, or exactly one arc with 1 or 0 supporting reads
                bool in = in_arc == lemon::INVALID || (arcs[in_arc] < 2 && ++in_arc == lemon::INVALID);
                bool out = out_arc == lemon::INVALID || (arcs[out_arc] < 2 && ++out_arc == lemon::INVALID);
                // a node is unique if for ingoing and outgoing arcs there is at most 1 supporting read
                if (in && out)
                    ++num_unique;
                else
                    ++num_non_unique;
            }
            viable = num_non_unique < 4 * num_unique;
        }
    }
    return *viable;
}

void DeBruijnGraph::prune(size_t threshold)
{
    // We need to store the nodes and arcs to be deleted, because the iterator gets confused otherwise.

    { // Erase each arc with read support below threshold.
        std::vector<lemon::ListDigraph::Arc> del;
        for (lemon::ListDigraph::ArcIt arc(graph); arc != lemon::INVALID; ++arc)
            if (arcs[arc] < threshold)
                del.push_back(arc);

        for (auto const & arc : del)
            graph.erase(arc);
    }

    { // Erase unconnected nodes.
        std::vector<lemon::ListDigraph::Node> del;
        for (lemon::ListDigraph::NodeIt node(graph); node != lemon::INVALID; ++node)
        {
            if (lemon::ListDigraph::OutArcIt(graph, node) == lemon::INVALID &&
                lemon::ListDigraph::InArcIt(graph, node) == lemon::INVALID)
                del.push_back(node);
        }

        for (auto const & node : del)
            graph.erase(node);
    }
}

std::vector<std::pair<float, seqan3::dna4_vector>> DeBruijnGraph::collect_haplotype_sequences() const
{
    // We do a depth-first search (dfs) in order to traverse all possible paths.
    lemon::Dfs<lemon::ListDigraph> dfs(graph);
    dfs.init();

    // Find all source nodes: they have no incoming arcs. Build a map of the out-degrees for each node.
    lemon::ListDigraph::NodeMap<size_t> out_degree(graph, 0);
    for (lemon::ListDigraph::NodeIt node(graph); node != lemon::INVALID; ++node)
    {
        if (lemon::ListDigraph::InArcIt(graph, node) == lemon::INVALID)
            dfs.addSource(node);
        for (lemon::ListDigraph::OutArcIt arc(graph, node); arc != lemon::INVALID; ++arc)
            out_degree[node] += arcs[arc];
    }

    // Storage for all the paths. Paths can be static, because only complete paths are stored here.
    std::vector<lemon::StaticPath<lemon::ListDigraph>> paths;
    {
        lemon::SimplePath<lemon::ListDigraph> path; // the currently constructed path
        lemon::ListDigraph::Node prev_node = graph.source(dfs.nextArc()); // the previously processed node

        while (!dfs.emptyQueue())
        {
            lemon::ListDigraph::Arc arc = dfs.processNextArc(); // go to next arc
            if (graph.source(arc) != prev_node) // we have jumped back to a fork
            {
                paths.emplace_back(path); // store the path that was completed before the jump
                while (!path.empty() && graph.target(path.back()) != graph.source(arc))
                    path.eraseBack(); // delete backwards until fork: the prefix is still needed for the next path
            }
            path.addBack(arc); // add to the current path
            prev_node = graph.target(arc); // update the processed node
        }
        paths.emplace_back(path); // store the last path that is completed when dfs finishes
        assert(lemon::checkPath(graph, path));
    }

    // For each path: collect sequence and score.
    std::vector<std::pair<float, seqan3::dna4_vector>> haplotypes(paths.size());
    for (auto && [path, haplotype] : seqan3::views::zip(paths, haplotypes))
    {
        auto && [score, seq] = haplotype;
        score = 1; // we will multiply the score, so start with 1
        seq.resize(path.length() + len); // allocate the haplotype sequence length

        // Convert the first node's kmer hash into sequence.
        size_t const hash = nodes[lemon::pathSource(graph, path)];
        for (unsigned char idx = 0; idx < len; ++idx)
            seq[idx].assign_rank(hash >> 2 * (len - idx - 1) & 0b11);

        // Iterate the path and collect further sequence characters.
        size_t idx = len;
        for (lemon::StaticPath<lemon::ListDigraph>::ArcIt arc(path); arc != lemon::INVALID; ++arc, ++idx)
        {
            seq[idx].assign_rank(nodes[graph.target(arc)] & 0b11);
            score *= static_cast<float>(arcs[arc]) / static_cast<float>(out_degree[graph.source(arc)]);
        }
    }
    std::sort(haplotypes.rbegin(), haplotypes.rend()); // highest score first
    return haplotypes;
}

void DeBruijnGraph::export_dot_format(std::string const & filename) const
{
    std::ofstream ofs(filename);
    ofs << "digraph G {\n";
    for (lemon::ListDigraph::NodeIt node(graph); node != lemon::INVALID; ++node)
    {
        // generate label
        std::string label;
        label.resize(len);
        for (unsigned char idx = 0; idx < len; ++idx)
            label[idx] = seqan3::dna4{}.assign_rank(nodes[node] >> 2 * (len - idx - 1) & 0b11).to_char();

        // produce output
        ofs << nodes[node] << " [label=\"" << label << "\"];\n";
        for (lemon::ListDigraph::OutArcIt arc(graph, node); arc != lemon::INVALID; ++arc)
            ofs << nodes[node] << " -> " << nodes[graph.runningNode(arc)] << " [label=\"" << arcs[arc] << "\"];\n";
    }
    ofs << "}\n";
}
