#pragma once

#include "Args.h"
#include "GeneNetwork.h"

#include <iosfwd>
#include <string>
#include <vector>
#include <unordered_map>

namespace gxna {

// Experiment is the main class that manages input data sets and runs the analysis.

template<typename Type> class MultipleTest;
class Permutation;

// Phenotype reads a sequence of sample labels and sequentially maps them to integers.
// For example, the sequence B B A C C A B A will be stored as 0 0 1 2 2 1 0 1,
// with nSamples() == 8 and nLabels() == 3.
class Phenotype {
 public:
    size_t nSamples() const { return m_sample.size(); }
    size_t nLabels() const { return m_label.size(); }
    const std::string name() const { return m_name; }
    const std::vector<int> get() const { return m_sample; }

    void read(std::istream&);
    friend std::ostream& operator<<(std::ostream&, const Phenotype&);
    void filter(const std::vector<std::string>& v);

 private:
    std::string m_name;
    std::vector<int> m_sample;
    std::vector<std::string> m_label;
    std::unordered_map<std::string, int> m_label2n;
};

struct GeneData {
    void printNameId(std::ostream& os) const;
    void print(std::ostream& os) const;
    std::string label() const;

    std::string id;  // NCBI gene ID
    std::string name = "NA";
    size_t nProbes = 0;  // probes that map to this gene AND have expression data
    std::vector<double> expression;
    double sd = 0;  // standard deviation of expression
    double shrinkageFactor = 1.0;  // used for empirical Bayes shrinkage
    double score = 0;
    double scorePerm = 0;  // score for the current permutation
};

class Experiment {
 public:
    explicit Experiment(Args&);
    void run();
    std::vector<double> operator()(const Permutation& perm);  // used by MultipleTest

 private:
    size_t nSamples() const { return m_phenotype[0].nSamples(); }

    // Computation
    double scoreNodeList(const GeneNetwork::NodeList& genes,
                         const std::vector<int>& pheno, int nLabels);
    void setShrinkageFactor();

    // Input
    size_t findPhenotype(const std::string& name) const;
    void readPhenotypes(const std::string& filename);
    void readProbes(const std::string& filename);
    void readGeneNames(const std::string& filename);
    void readExpression(const std::string& filename);

    // Output
    void writeHTML(const std::string& htmlFilename, const std::string& frameFilename,
                   const std::string& startingFrame) const;
    void printResults(const MultipleTest<Experiment>&);

    struct TestData {
        size_t root;
        GeneNetwork::NodeList cluster;
        double score = 0;
    };

    Args& args;  // program arguments
    std::vector<GeneData> m_gene;
    std::unordered_map<std::string, std::string> m_probe2geneID;
    std::unordered_map<std::string, size_t> m_geneID2index;
    GeneNetwork m_geneNetwork;  // gene interaction graph
    std::vector<Phenotype> m_phenotype;
    Phenotype m_mainPhenotype;  // phenotype used to compute scores
    std::vector<TestData> m_testData;  // one for each root/hypothesis being tested
    size_t m_permCount = 0;
};

}  // namespace gxna
