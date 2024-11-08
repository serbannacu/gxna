#pragma once

#include <iosfwd>
#include <string>
#include <vector>

namespace gxna {

// Program parameters

enum class AlgoType { Basic, GXNA };

std::ostream& operator<<(std::ostream& os, const AlgoType& x);

struct Args {
    Args();
    void check();
    void parse(int argc, char *argv[]);  // read from command line
    void read(const std::string& filename, bool strict);  // read from file
    void print(std::ostream&) const;
    void usage(std::ostream&) const;  // help message

    // Filenames
    std::string name;
    std::string version;
    std::string refDir;  // reference data
    std::string expDir;  // experiment data (expression, phenotypes)
    std::string outputDir;
    std::string geneFile;  // gene metadata
    std::string interactionFile;  // gene interactions
    std::string probeFile;  // probe metadata
    std::string expressionFile;
    std::string phenotypeFile;

    // Algorithm phenotypes
    std::string test;  // test this phenotype
    std::string invariant;  // use this phenotype for invariant permutations
    std::vector<std::string> filter;  // filter for these phenotype labels

    // Algorithm parameters
    AlgoType algoType;
    int radius;  // ball size for Basic algo
    int depth;  // subgraph size for GXNA algo
    bool flexSize;
    double minSD;
    int minDegree;
    int minDistance;

    // Algorithm subgraph score calculations
    bool sumScore;
    bool sumSigned;
    double scalingExponent;

    // Algorithm p-value calculations
    bool maxTscaled;
    int nPerms;  // for resampling

    // Algorithm other parameters
    bool shrink;  // adjust gene scores using empirical Bayes shrinkage
    int seed;  // for random number generator

    // Output
    int nRows;  // max number of graphs to display in html file
    int nDetailed;  // max number of graphs to write or draw details
    double maxOverlap;
    bool draw;  // use Graphviz to render graphs
    bool progress;  // show progress indicator

 private:
    void setFilenames();
    bool setImpl(const std::string& key, const std::string& val);
    void set(const std::string& key, const std::string& val);
};

}  // namespace gxna

