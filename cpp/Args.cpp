#include "Args.h"
#include "VectorUtil.h"
#include "Exception.h"

#include <fstream>
#include <iostream>
#include <iomanip>
#include <sstream>
#include <stdexcept>

namespace gxna {

std::ostream& operator<<(std::ostream& os, const AlgoType& x) {
    switch (x) {
    case AlgoType::Basic:
        os << "Basic";
        break;
    case AlgoType::GXNA:
        os << "GXNA";
        break;
    default:
        throw Exception("Bad AlgoType ") << int(x);
    }
    return os;
}

Args::Args()
    : version("000"),
      refDir("refdata"),
      expDir("expdata"),
      outputDir("output"),
      geneFile("geneID2name"),
      interactionFile("human.gra"),
      algoType(AlgoType::Basic),
      radius(0),
      depth(15),
      flexSize(false),
      minSD(0),
      minDegree(0),
      minDistance(0),
      sumScore(true),
      sumSigned(true),
      scalingExponent(0.6),
      maxTscaled(false),
      nPerms(100),
      shrink(false),
      seed(5),
      nRows(250),
      nDetailed(25),
      maxOverlap(0.75),
      draw(false),
      progress(true)
{}

void Args::check() {
    // Filter out genes without interactions, unless doing single gene analysis.
    if (minDegree == 0 && (algoType == AlgoType::GXNA || radius > 0))
        minDegree = 1;

    const char *msg = nullptr;
    if (!name.size())
        msg = "Empty name";
    else if (!probeFile.size())
        msg = "Must set probeFile";
    else if (!sumScore && algoType == AlgoType::GXNA)
        msg = "Need sumScore = true for algoType = GXNA";
    else if (!sumScore && shrink)
        msg = "Need sumScore = true to use shrinkage";
    else if (filter.size() == 1)
        msg = "Filter must be empty or have at least two values";
    if (msg)
        throw Exception(msg);
}

void Args::parse(int argc, char *argv[]) {
    if (argc % 2 == 0)
        throw Exception("Each key token must be followed by a value token");
    for (int n = 1; n <= argc - 2; n += 2) {
        char *s = argv[n];
        if (s[0] == '-') {
            std::string key(++s);
            std::string val(argv[n+1]);
            if (key == "argFile")
                read(val, true);  // exit if file is missing
            else
                set(key, val);
        }
        else {
            throw Exception(std::string("Bad command line token ") + s);
        }
    }
}

void Args::read(const std::string& filename, bool strict) {
    std::ifstream is(filename.c_str());
    if (is)
        std::cout << "Reading args from " << filename << std::endl;
    else if (strict)
        throw Exception("Could not open " + filename);
    std::string key, val;
    while (is >> key >> val)
        set(key, val);
}

void Args::print(std::ostream& os) const {
    os << std::boolalpha;
    os << "name " << name << '\n';
    os << "version " << version << '\n';
    os << "refDir " << refDir << '\n';
    os << "expDir " << expDir << '\n';
    os << "outputDir " << outputDir << '\n';
    os << "geneFile " << geneFile << '\n';
    os << "interactionFile " << interactionFile << '\n';
    os << "probeFile " << probeFile << '\n';
    os << "expressionFile " << expressionFile << '\n';
    os << "phenotypeFile " << phenotypeFile << '\n';
    os << "test " << test << '\n';
    os << "invariant " << invariant << '\n';
    os << "filter " << filter << '\n';
    os << "algoType " << algoType << '\n';
    os << "radius " << radius << '\n';
    os << "depth " << depth << '\n';
    os << "flexSize " << flexSize << '\n';
    os << "minSD " << minSD << '\n';
    os << "minDegree " << minDegree << '\n';
    os << "minDistance " << minDistance << '\n';
    os << "sumScore " << sumScore << '\n';
    os << "sumSigned " << sumSigned << '\n';
    os << "scalingExponent " << scalingExponent << '\n';
    os << "maxTscaled " << maxTscaled << '\n';
    os << "nPerms " << nPerms << '\n';
    os << "shrink " << shrink << '\n';
    os << "seed " << seed << '\n';
    os << "nRows " << nRows << '\n';
    os << "nDetailed " << nDetailed << '\n';
    os << "maxOverlap " << maxOverlap << '\n';
    os << "draw " << draw << '\n';
    os << "progress " << progress << '\n';
}

static void log_option(std::ostream& os, const char *name, const char *arg,
                       const char *description) {
    os << "  -" << std::left
       << std::setw(16) << name << ' '
       << std::setw(8) << arg << "  "
       << std::right << description << '\n';
}

void Args::usage(std::ostream& os) const {
    os << "Usage: gxna [options]\n";
    os << "\n";
    os << "Options:\n";
    log_option(os, "name", "string", "experiment name");
    log_option(os, "version", "string", "output version");
    log_option(os, "interactionFile", "filename", "gene interaction file");
    log_option(os, "probeFile", "filename", "probe annotation file");
    log_option(os, "test", "string", "test this phenotype");
    log_option(os, "invariant", "string", "use this phenotype for invariant permutations");
    log_option(os, "algoType", "type", "Basic | GXNA");
    log_option(os, "radius", "int", "ball radius for ball search");
    log_option(os, "depth", "int", "max cluster size for adapted search");
    log_option(os, "nPerms", "int", "number of permutations");
    log_option(os, "shrink", "bool", "adjust scores via empirical Bayes shrinkage");
    log_option(os, "seed", "int", "seed for random number generator");
    log_option(os, "nRows", "int", "max number of graphs to display in html file");
    log_option(os, "nDetailed", "int", "max number of graphs to write or draw details");
    log_option(os, "draw", "bool", "render graphs to SVG files (needs Graphviz)");
    os << "\n";
}

void from_string(std::string& lhs, const std::string& rhs) {
    lhs = rhs;
}

void from_string(double& lhs, const std::string& rhs) {
    lhs = std::stod(rhs);
}

void from_string(int& lhs, const std::string& rhs) {
    lhs = std::stoi(rhs);
}

void from_string(bool& lhs, const std::string& rhs) {
    if (rhs == "1" || rhs == "T" || rhs == "true" || rhs == "True" || rhs == "TRUE")
        lhs = true;
    else if (rhs == "0" || rhs == "F" || rhs == "false" || rhs == "False" || rhs == "FALSE")
        lhs = false;
    else
        throw Exception("Bad bool value " + rhs);
}

template<typename T>
void from_string(std::vector<T>& lhs, const std::string& rhs) {
    std::istringstream is(rhs);
    is >> lhs;
}

void from_string(AlgoType& lhs, const std::string& rhs) {
    if (rhs == "basic" || rhs == "Basic" || rhs == "BASIC")
        lhs = AlgoType::Basic;
    else if (rhs == "gxna" || rhs == "Gxna" || rhs == "GXNA")
        lhs = AlgoType::GXNA;
    else
        throw Exception("Bad AlgoType " + rhs);
}

void Args::setFilenames() {
    if (!expressionFile.size())
        expressionFile = name + ".exp";
    if (!phenotypeFile.size())
        phenotypeFile = name + ".phe";
    auto argsFile = name + ".arg";
    read(expDir + "/" + argsFile, false);  // ignore if file is missing
}

bool Args::setImpl(const std::string& key, const std::string& val) {
    if (key == "name") {
        from_string(name, val);
        setFilenames();
    }
    else if (key == "version")
        from_string(version, val);
    else if (key == "refDir")
        from_string(refDir, val);
    else if (key == "expDir")
        from_string(expDir, val);
    else if (key == "outputDir")
        from_string(outputDir, val);
    else if (key == "geneFile")
        from_string(geneFile, val);
    else if (key == "interactionFile")
        from_string(interactionFile, val);
    else if (key == "probeFile")
        from_string(probeFile, val);
    else if (key == "expressionFile")
        from_string(expressionFile, val);
    else if (key == "phenotypeFile")
        from_string(phenotypeFile, val);
    else if (key == "test")
        from_string(test, val);
    else if (key == "invariant")
        from_string(invariant, val);
    else if (key == "filter")
        from_string(filter, val);
    else if (key == "algoType")
        from_string(algoType, val);
    else if (key == "radius")
        from_string(radius, val);
    else if (key == "depth")
        from_string(depth, val);
    else if (key == "flexSize")
        from_string(flexSize, val);
    else if (key == "minSD")
        from_string(minSD, val);
    else if (key == "minDegree")
        from_string(minDegree, val);
    else if (key == "minDistance")
        from_string(minDistance, val);
    else if (key == "sumScore")
        from_string(sumScore, val);
    else if (key == "sumSigned")
        from_string(sumSigned, val);
    else if (key == "scalingExponent")
        from_string(scalingExponent, val);
    else if (key == "maxTscaled")
        from_string(maxTscaled, val);
    else if (key == "nPerms")
        from_string(nPerms, val);
    else if (key == "shrink")
        from_string(shrink, val);
    else if (key == "seed")
        from_string(seed, val);
    else if (key == "nRows")
        from_string(nRows, val);
    else if (key == "nDetailed")
        from_string(nDetailed, val);
    else if (key == "maxOverlap")
        from_string(maxOverlap, val);
    else if (key == "draw")
        from_string(draw, val);
    else if (key == "progress")
        from_string(progress, val);
    else
        return false;
    return true;
}

void Args::set(const std::string& key, const std::string& val) {
    try {
        if (!setImpl(key, val))
            throw Exception("Unknown key " + key);
    }
    catch (const std::invalid_argument& ex) {
        throw Exception("Bad key/val " + key + " = " + val);
    }
}

}  // namespace gxna
