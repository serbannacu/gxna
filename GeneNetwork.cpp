#include "GeneNetwork.h"
#include "Exception.h"
#include "Statistics.h"
#include "VectorUtil.h"

#include <algorithm>
#include <fstream>
#include <iostream>
#include <sstream>
#include <unordered_set>

namespace gxna {

// Read gene network from file
// Each line represents a gene interaction, in the format: gene1 gene2 [type [source]]
// A gene pair could have multiple interaction types/sources but only the first one is processed

void GeneNetwork::readInteractions(const std::string& filename, std::unordered_map<std::string, int>& gene2index) {
    std::ifstream is(filename.c_str());
    if (!is)
        throw Exception("Could not open " + filename);
    std::string line;
    while (getline(is, line)) {
        std::istringstream ss(line);
        std::string gene1, gene2, type, source;
        ss >> gene1 >> gene2;
        if (!ss) {
            std::cerr << "readInteractions bad line " << line << '\n';
        }
        else {
            auto p1 = gene2index.find(gene1);
            auto p2 = gene2index.find(gene2);
            if (p1 != gene2index.end() && p2 != gene2index.end()) { // ignore unknown genes
                int n1 = p1->second;
                int n2 = p2->second;
                ss >> type >> source; // not using source for now
                if (n1 != n2) // ignore self-interactions
                    addEdge(n1, n2, type);
            }
        }
    }
    std::cerr << "Read " << m_nEdges << " interactions from " << filename << '\n';
    for (auto it : gene2index)
        setLabel(it.second, it.first);
}

void GeneNetwork::addEdge(int n1, int n2, const std::string& type) {
    Edge edge(n1, n2);
    if (m_edgeType.find(edge) == m_edgeType.end()) { // new edge
        if (m_edgeType.find(Edge(n2, n1)) == m_edgeType.end()) {
            auto n = std::max(n1, n2);
            if (m_neighbors.size() <= n)
                m_neighbors.resize(n + 1);
            m_neighbors[n1].push_back(n2);
            m_neighbors[n2].push_back(n1);
            m_nEdges++;
        }
        m_edgeType[edge] = type;
    }
}

GeneNetwork::NodeList GeneNetwork::ball(int center, int radius) const {
    std::unordered_set<int> ball, boundary;
    ball.insert(center);
    boundary.insert(center);
    for (int depth = 1; depth <= radius; ++depth) {
        std::unordered_set<int> newBoundary;
        for (int v : boundary)
            for (int w : m_neighbors[v])
                if (ball.find(w) == ball.end())
                    newBoundary.insert(w);
        boundary = newBoundary;
        ball.insert(boundary.begin(), boundary.end());
    }
    return NodeList(ball.begin(), ball.end());
}

void GeneNetwork::setScores(const std::vector<double>& score, double scalingExponent) {
    assert(score.size() == nNodes());
    FastDataSet ds(score);
    m_score = score;
    m_taken.resize(nNodes());
    m_meanScore = ds.mean();
    m_scalingExponent = scalingExponent;
    for (auto& nodes : m_neighbors) // sort neighbors by decreasing score
        std::sort(nodes.begin(), nodes.end(), [&](int j, int k) { return m_score[j] > m_score[k]; });
}

double GeneNetwork::subgraphScore(const NodeList& subgraph) const {
    double sum = 0;
    for (int node : subgraph)
        sum += m_score[node];
    return getScaledScore(sum, subgraph.size());
}

// Main search algorithm: finds high scoring connected subgraph by greedy expansion starting at root
// Before using it, must call setScores, which pre-sorts each node's neighbors
// At each step, findSubgraph does a linear search over the current subgraph and adds the neighbor with the highest score
// In our main use case (sparse graph, small subgraphs) this performs better empirically than a priority queue

double GeneNetwork::findSubgraph(int root, int depth, bool flexSize, NodeList& result) {
    std::vector<int> subgraph;
    std::vector<int> index(depth); // for each node, an index to its best scoring neighbor not yet added to the subgraph
    subgraph.reserve(depth);
    subgraph.emplace_back(root);
    m_taken[root] = true;
    double sumScore = m_score[root];

    while (subgraph.size() < depth) {
        double bestScore = 0; // only add nodes with positive scores
        int node, bestNode = -1;
        for (int i = 0; i < subgraph.size(); ++i) { // for each current node find its best neighbor
            auto& neighbors = m_neighbors[subgraph[i]]; // pre-sorted with the highest score at the front
            int j = index[i];
            while (j < neighbors.size() && m_taken[node = neighbors[j]])
                ++j;
            index[i] = j;
            if (j < neighbors.size()) {
                double score = m_score[node];
                if (bestScore < score) {
                    bestScore = score;
                    bestNode = node;
                }
            }
        }
        if (bestNode == -1) // there are no more nodes to be added
            break;
        else {
            subgraph.emplace_back(bestNode);
            m_taken[bestNode] = true;
            sumScore += bestScore;
        }
    }
    for (int node : subgraph) // clean up m_token
        m_taken[node] = false;

    double scaledScore = getScaledScore(sumScore, subgraph.size());
    if (flexSize) { // check if some sub-subgraph has a better (scaled) score than the subgraph
        int scaledSize = subgraph.size();
        double val = 0;
        for (int size = 1; size < subgraph.size(); ++size) {
            val += m_score[subgraph[size - 1]];
            double score = getScaledScore(val, size);
            if (scaledScore < score) {
                scaledScore = score;
                scaledSize = size;
            }
        }
        if (scaledSize < subgraph.size())
            subgraph.resize(scaledSize);
    }
    result = subgraph;
    return scaledScore;
}

void GeneNetwork::print(std::ostream& os) const {
    os << nNodes() << " nodes " << m_nEdges << " edges\n";
    for (int v = 0; v < nNodes(); ++v)
        if (degree(v))
            os << v << ' ' << degree(v) << " : " << m_neighbors[v] << '\n';
}

void GeneNetwork::write(const NodeList& subgraph, const std::string& filenameDOT, const std::string& filenameSVG) const {
    std::ofstream osDOT(filenameDOT.c_str());
    writeDOT(subgraph, osDOT);
    osDOT.flush();
    if (filenameSVG.size()) {
        std::string cmd = "neato -Tsvg -o" + filenameSVG + " " + filenameDOT;
        auto val = system(cmd.c_str());
        if (val)
            std::cerr << "Command " << cmd << " failed with exit value " << val << '\n';
    }
}

static void writeDOTNode(std::ostream& os, int v, const std::string& label, const std::string& color = "") {
    os << v << " [label = \"" << label << "\"";
    if (color.size())
        os << ", color=\"" + color + "\"";
    os << " ] ; \n";
}

void GeneNetwork::writeDOTEdge(std::ostream& os, int v, int w) const {
    auto it = m_edgeType.find(Edge(v, w));
    if (it != m_edgeType.end()) {
        auto& type = it->second;
        os << v << " -> " << w;
        if (type == "I") // undirected interaction
            os << " [arrowhead = \"none\" ]";
        else
            os << " [label = \"" << type << "\" ]";
        os << '\n';
    }
}

void GeneNetwork::writeDOT(const NodeList& subgraph, std::ostream& os) const {
    os << "digraph G {\n"; // print header
    os << "overlap = scale ;\n";
    for (auto v : subgraph)
        writeDOTNode(os, v, v < m_label.size() ? m_label[v] : "");
    std::unordered_set<int> nodes(subgraph.begin(), subgraph.end());
    for (auto v : subgraph) {
        for (auto w : m_neighbors[v]) {
            if (v < w && nodes.find(w) != nodes.end()) {
                writeDOTEdge(os, v, w);
                writeDOTEdge(os, w, v);
            }
        }
    }
    os << "}\n"; // print trailer
}

} // namespace gxna
