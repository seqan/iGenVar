#pragma once

#include <optional>
#include <string>

#include <seqan3/alphabet/nucleotide/dna4.hpp>

#include <lemon/maps.h>
#include <lemon/list_graph.h>

//! \brief A De Bruijn Graph represents overlaps between DNA sequences.
class DeBruijnGraph
{
private:
    unsigned char len; //!> The length of the sequence fragments in the nodes (k-mer length).
    std::optional<bool> viable; //!> Flag whether the graph is complex and acyclic. No value means unknown.
    lemon::ListDigraph graph; //!> The graph structure.
    lemon::ListDigraph::ArcMap<size_t> arcs; //!> Each arc is assigned the number of reads that support it.
    lemon::CrossRefMap<lemon::ListDigraph, lemon::ListDigraph::Node, size_t> nodes; //!> Each node represents a k-mer.

public:
    /*!\name Constructor and destructor
     * \{
     */
    DeBruijnGraph() : len(0), viable(false), graph(), arcs(graph), nodes(graph) {} //!< Defaulted.
    DeBruijnGraph(DeBruijnGraph const &) = delete; //!< Deleted.
    DeBruijnGraph(DeBruijnGraph &&) noexcept = delete; //!< Deleted.
    DeBruijnGraph & operator=(DeBruijnGraph const &) = delete; //!< Deleted.
    DeBruijnGraph & operator=(DeBruijnGraph &&) noexcept = delete; //!< Deleted.
    ~DeBruijnGraph() = default; //!< Defaulted.
    //!\}

    /*!
     * \brief Reset and initialize the graph with a reference sequence.
     * \param[in] kmer_len The k-mer length.
     * \param[in] reference The reference sequence that is added to the graph.
     * \return whether the reference k-mers are unique. If false, please re-initialize with higher `len`.
     */
    bool init_sequence(unsigned char kmer_len, seqan3::dna4_vector const & reference);

    /*!
     * \brief Add a read sequence to the graph. The `viable` property will be unknown afterwards.
     * \param[in] read The read sequence.
     */
    void add_read(seqan3::dna4_vector const & read);

    /*!
     * \brief Check if the graph is acyclic (DAG) and complex.
     * \return whether the graph is viable.
     *
     * \details
     * Complex means: number of non-unique nodes < 4 * number of unique nodes.
     * If the `viable` property is known, it is returned, otherwise it is computed once.
     */
    bool is_viable();

    /*!
     * \brief Prune arcs that are not well supported by reads. Afterwards, remove nodes that do not have an arc.
     * \param[in] threshold Arcs with fewer than `threshold` supporting reads are removed.
     */
    void prune(size_t threshold);

    /*!
     * \brief Traverse the graph and collect the haplotype sequences.
     * \return a vector of haplotype sequences with associated score, sorted by score in descending order.
     *
     * \details
     * The scores add up to 1 and represent the relative read support of the path.
     */
    std::vector<std::pair<float, seqan3::dna4_vector>> collect_haplotype_sequences() const;

    /*!
     * \brief Write the graph in dot format to a file. Can be used to visualize the graph with the dot program.
     * \param[in] filename The filename of the graph file (usually with suffix .gv).
     */
    void export_dot_format(std::string const & filename) const;
};
