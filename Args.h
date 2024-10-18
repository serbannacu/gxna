#pragma once

#include <iosfwd>
#include <string>
#include <vector>

namespace gxna {

enum class AlgoType { Basic, GXNA }; 

std::ostream& operator<<(std::ostream& os, const AlgoType& x);

struct Args {
    Args();
    void check();
    void parse(int argc, char *argv[]);
    void read(const std::string& filename, bool strict = true);
    void print(std::ostream&) const;

    // filenames
    std::string name;
    std::string version;
    std::string refDir;
    std::string inputDir;
    std::string outputDir;
    std::string geneFile;
    std::string interactionFile;
    std::string probeFile;
    std::string expressionFile;
    std::string phenotypeFile;
    std::string typeFile;

    // algorithm filter and search
    std::vector<std::string> phenotypes; // only use these phenotypes
    AlgoType algoType;
    int radius; // ball size for Basic algo, also used for filtering 
    int depth; // subgraph size for GXNA algo
    bool flexSize;
    double minSD;
    int minDegree;

    // algorithm subgraph score calculations
    bool sumScore;
    bool sumSigned;
    double scalingExponent;

    // algorithm p-value calculations
    bool maxTscaled;
    int nPerms;
    bool useInvariantPerms;

    // algorithm other parameters
    bool shrink; // adjust gene scores using empirical Bayes shrinkage
    int seed; // for random number generator

    // output
    int graphCount; // max number of graphs to draw/write details
    double maxOverlap;
    int maxRows; // max number of graphs to display in html file
    bool draw; // use Graphviz to render graphs

private:
    void setFilenames();
    bool setImpl(const std::string& key, const std::string& val);
    void set(const std::string& key, const std::string& val);
};

} // namespace gxna

