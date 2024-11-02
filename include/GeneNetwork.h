#pragma once

#include <cmath>
#include <iosfwd>
#include <string>
#include <vector>
#include <unordered_map>

// A gene network is an undirected (sparse) graph where
// nodes represent genes and edges represent interactions.
// Nodes are indexed by integers, and each node stores the list of its neighbors.
// Interaction type data is stored as a string but only used for drawing the network.

namespace gxna {

class GeneNetwork {
 public:
    typedef std::pair<size_t, size_t> Edge;
    typedef std::vector<size_t> NodeList;

    void readInteractions(const std::string& filename,
                          std::unordered_map<std::string, size_t>& gene2index);

    void setSize(size_t n) {
        m_neighbors.resize(n);
    }

    void setLabel(size_t n, const std::string& label) {
        if (n >= m_label.size())
            m_label.resize(n + 1);
        m_label[n] = label;
    }

    size_t nNodes() const {
        return m_neighbors.size();
    }

    size_t degree(size_t node) const {
        return m_neighbors[node].size();
    }

    // Compute list of nodes at distance <= radius from root
    NodeList ball(size_t root, size_t radius) const;

    // Set nodes scores: must be used before calling findSubgraph
    void setScores(const std::vector<double>& score, double scalingExponent);

    double subgraphScore(const NodeList& subgraph) const;

    // Find highest scoring subgraph for a given root, using greedy search
    double findSubgraph(size_t root, size_t depth, bool flexSize, NodeList& result);

    // Print all nodes and edges
    void print(std::ostream& os) const;

    // Write subgraph into Graphviz DOT file, and optionally render into SVG file
    void write(const NodeList& subgraph, const std::string& filenameDOT,
               const std::string& filenameSVG = "") const;

 private:
    void addEdge(size_t n1, size_t n2, const std::string& type);

    // Compute scaled score of a subgraph, given the sum of its node scores and size
    double getScaledScore(double sumScore, size_t size) const {
        return (sumScore - size * m_meanScore) / std::pow(size, m_scalingExponent);
    }

    void writeDOTEdge(std::ostream& os, size_t v, size_t w) const;
    void writeDOT(const NodeList& subgraph, std::ostream& os) const;

    // needed for unordered_map
    struct pair_hash {
        template <class T1, class T2>
        std::size_t operator () (const std::pair<T1, T2> &p) const {
            auto x1 = std::hash<T1>{}(p.first);
            auto x2 = std::hash<T2>{}(p.second);
            x1 ^= x2 + 0x9e3779b9 + (x1 << 6) + (x1 >> 2);  // boost::hash_combine
            return x1;
        }
    };

    std::vector<NodeList> m_neighbors;  // for each vertex a list of neighbors
    std::vector<std::string> m_label;  // (optional) node labels, used for output
    std::unordered_map<Edge, std::string, pair_hash> m_edgeType;  // used for output
    size_t m_nEdges = 0;

    std::vector<double> m_score;  // a score for each node
    std::vector<bool> m_taken;  // a flag for each node, true if node is in current subgraph
    double m_meanScore = 0;
    double m_scalingExponent = 0;
};

}  // namespace gxna
