// Master class that manages the various data sets and runs the analysis

#pragma once

#include "Args.h"
#include "GeneNetwork.h"

#include <iosfwd>
#include <string>
#include <vector>
#include <unordered_map>

namespace gxna {

template<typename Type> class MultipleTest;
class Permutation;

class PhenotypeVector {
public:
    size_t nSamples() const { return m_type.size(); }
    size_t nTypes() const { return m_type2name.size(); }
    const std::vector<int>& type() const { return m_type; }
    
    void read(const std::string& filename);
    void filter(const std::vector<std::string>& v);

private:
    std::vector<int> m_type;
    std::vector<std::string> m_type2name;
    std::unordered_map<std::string, int> m_name2type;
};

struct GeneData {
    void print(std::ostream& os) const;
    std::string label() const;

    std::string id; // NCBI gene ID
    std::string name = "NA";
    int nProbes = 0; // number of probes that map to this gene AND have expression data
    std::vector<double> expression;
    double sd = 0; // standard deviation of expression
    double shrinkageFactor = 1.0; // used for moderated test statistics (empirical Bayes shrinkage)
    double score = 0;
    double scorePerm = 0; // score for the current permutation
};

class Experiment {
public:
    Experiment(Args&);
    void run();
    std::vector<double> operator()(const Permutation& perm); // used by MultipleTest

private:
    // Computation
    void setShrinkageFactor();
    double scoreNodeList(const GeneNetwork::NodeList& genes, const std::vector<int>& phenotype);
    
    // Input
    void readProbes(const std::string& filename);
    void readGeneNames(const std::string& filename);
    void readExpression(const std::string& filename);

    // Output
    void writeHTML(const std::string& htmlFilename, const std::string& frameFilename, const std::string& startingFrame) const;
    void printResults(const MultipleTest<Experiment>&);

    struct TestData {
        int root;
        GeneNetwork::NodeList cluster;
        double score = 0;
    };
    
    Args& args; // program arguments
    std::vector<GeneData> m_gene; 
    std::unordered_map<std::string, std::string> m_probe2geneID;
    std::unordered_map<std::string,int> m_geneID2index;
    GeneNetwork m_geneNetwork; // gene interaction graph
    PhenotypeVector m_phenotype; // main phenotype
    PhenotypeVector m_phenotypeInvariant; // only used to generate invariant permutations
    std::vector<TestData> m_testData; // one for each root/hypothesis being tested
    int m_permCount = 0;
};

} // namespace gxna
