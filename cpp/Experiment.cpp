#include "Experiment.h"

#include "Distribution.h"
#include "Exception.h"
#include "Log.h"
#include "MultipleTest.h"
#include "RandomPermutation.h"
#include "Statistics.h"
#include "VectorUtil.h"

#include <filesystem>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <unordered_set>

namespace gxna {

// Read one phenotype from stream.
// Format: name label_1 ... label_m (one label per sample)
void Phenotype::read(std::istream& is) {
    is >> m_name;
    if (!is)
        throw Exception("Phenotype needs name");
    std::string label;
    while (is >> label) {
        if (label == m_name)
            throw Exception("Phenotype") << ' ' << m_name << ": sample label same as name";
        auto it = m_label2n.find(label);
        int i;
        if (it == m_label2n.end()) {  // new label
            m_label2n[label] = i = m_label.size();
            m_label.emplace_back(label);
        }
        else
            i = it->second;
        m_sample.emplace_back(i);
    }
    if (m_sample.empty())
        throw Exception("Phenotype") << ' ' << m_name << ": needs at least one sample";
    if (nLabels() < 2)
        throw Exception("Phenotype") << ' ' << m_name << ": needs at least two labels";
}

// Print name, labels and label counts.
std::ostream& operator<<(std::ostream& os, const Phenotype& p) {
    std::vector<int> count(p.nLabels());
    for (auto& val : p.m_sample)
        ++count[val];
    os << p.name() << " {";
    for (size_t i = 0; i < count.size(); ++i) {
        if (i)
            os << ", ";
        os << p.m_label[i] << " : " << count[i];
    }
    os << '}';
    return os;
}

void Phenotype::filter(const std::vector<std::string>& v) {
    for (auto& label : v) {
        auto it = m_label2n.find(label);
        if (it == m_label2n.end())
            throw Exception("Phenotype") << ' ' << m_name << ": unknown label " << label;
        // should use it->second
    }
    throw Exception("Phenotype::filter not implemented");
}

void GeneData::printNameId(std::ostream& os) const {
    os << std::setw(10) << std::left << name << ' '
       << std::setw(6) << std::right << id << ' ';
}

void GeneData::print(std::ostream& os) const {
    printNameId(os);
    os << std::fixed << std::setprecision(3)
       << std::setw(2) << nProbes << ' '
       << std::setw(7) << sd << ' '
       << std::setw(7) << score;
}

static std::ostream& operator<<(std::ostream& os, const GeneData& x) {
    x.print(os);
    return os;
}

std::string GeneData::text() const {
    std::ostringstream os;
    os << name << std::fixed << std::setprecision(2) << "\\n" << score;
    return os.str();
}

Experiment::Experiment(Args& args_)
    : args(args_) {

    readPhenotypes(args.expDir + "/" + args.phenotypeFile);
    readProbes(args.refDir + "/" + args.probeFile);
    readGeneNames(args.refDir + "/" + args.geneFile);
    readExpression(args.expDir + "/" + args.expressionFile);
    m_geneNetwork.setSize(m_gene.size());  // so genes with no interactions are included
    m_geneNetwork.readInteractions(args.refDir + "/" + args.interactionFile, m_geneID2index);
}

void Experiment::run() {
    // Prepare output directory and log stream.

    std::string outputPath = args.outputDir + "/" + args.name + "/" + args.version;
    std::filesystem::create_directories(outputPath);
    LogGuard logGuard {(outputPath + "/" + "log.txt").c_str()};  // use clog for output

    // Compute gene statistics.

    for (auto& data : m_gene) {
        FastDataSet d(data.expression);
        data.sd = d.sd();
    }
    if (args.shrink)  // use empirical Bayes shrinkage
        setShrinkageFactors();

    // Select root genes for multiple testing.
    // Potential roots will not be within args.minDistance of each other.

    std::vector<bool> ignore(m_gene.size());
    for (size_t root = 0; root < m_geneNetwork.nNodes(); ++root) {
        if (m_geneNetwork.degree(root) >= unsigned(args.minDegree)
            && m_gene[root].sd >= args.minSD && !ignore[root]) {  // valid root
            Cluster cluster;
            cluster.root = root;
            if (args.algoType == AlgoType::Basic)  // cluster is preset
                cluster.nodes = m_geneNetwork.ball(root, args.radius);
            if (args.minDistance > 0) {  // do not use nearby genes as roots
                for (auto v : m_geneNetwork.ball(root, args.minDistance))
                    ignore[v] = true;
            }
            m_cluster.emplace_back(cluster);
        }
    }

    // Prepare permutations for resampling.

    Permutation::seed(args.seed);
    PermutationGenerator *pg;
    if (!args.invariant.empty()) {
        auto& phenotype = m_phenotype[findPhenotype(args.invariant)];
        std::cout << "Using invariant permutations for phenotype "
                  << phenotype << std::endl;
        pg = new InvariantPermutation(nSamples(), args.nPerms, phenotype.get());
    }
    else
        pg = new UniformPermutation(nSamples(), args.nPerms);
    pg->setVerbose(args.progress);

    // Multiple roots can yield the same subgraph but with different node order.
    // For example, if x, y and z are genes connected only to each other,
    // they may generate the subgraphs {x, y, z}, {y, x, z} and {z, x, y}.
    // The subgraphs should have the same score, but they may not, since
    // floating point sums are not exactly associative due to rounding errors.

    // This can cause inconsistent behavior across platforms.
    // To prevent this, we provide a tolerance parameter to MultipleTest,
    // small but expected to exceed the rounding error.

    constexpr double Eps = 1e-9;
    MultipleTest<Experiment> mt(*this, m_cluster.size(), Eps);
    std::cout << "Testing " << m_cluster.size() << " objects using phenotype "
              << m_mainPhenotype << std::endl;
    mt.test(*pg, args.maxTscaled);
    delete pg;
    printResults(mt, outputPath);
}

static double computeScore(const std::vector<double>& expression,
                           const std::vector<int>& label, int nLabels) {
    double score;
    if (nLabels > 2) {  // F statistic
        double f = fstatLabel(&expression[0], label, nLabels);
        score = f2z(f, nLabels - 1, label.size() - nLabels);
    }
    else {  // t statistic; may want to convert to z-score
        double t = tstatLabel(&expression[0], label, 0, 1);  // pheno 0 vs 1
        score = t;
    }

    // The F and t statistics could in principle be infinite.
    // This would cause problems when scaling and/or comparing scores.
    // In real data, extremely large values are likely overly optimistic anyway.
    // Therefore, we enforce a bound, larger than any realistic z-score.

    constexpr double MaxScore = 35;
    if (std::fabs(score) > MaxScore)
        score = score < 0 ? -MaxScore : MaxScore;
    return score;
}

std::vector<double> Experiment::operator()(const Permutation& perm) {
    // Recompute gene scores.

    auto label = perm.apply(m_mainPhenotype.get());
    auto nLabels = m_mainPhenotype.nLabels();

    std::vector<double> geneScorePermAbs;  // one score per gene
    for (auto& data : m_gene) {
        double score = computeScore(data.expression, label, nLabels);
        if (args.shrink)
            score *= data.shrinkageFactor;
        data.scorePerm = score;
        if (!m_permCount)  // first perm is Id so labels are not scrambled
            data.score = score;  // save true score
        geneScorePermAbs.emplace_back(std::fabs(score));
    }
    m_geneNetwork.setScores(geneScorePermAbs, args.scalingExponent);

    // Recompute cluster scores.

    std::vector<double> clusterScorePermAbs;  // one score per cluster
    clusterScorePermAbs.reserve(m_cluster.size());
    if (args.algoType == AlgoType::Basic) {
        for (auto& cluster : m_cluster) {
            auto score = scoreNodeList(cluster.nodes, label, nLabels);
            clusterScorePermAbs.emplace_back(std::fabs(score));
            if (!m_permCount)
                cluster.score = score;
        }
    }
    else {  // AlgoType::GXNA
        GeneNetwork::NodeList nodes;
        for (auto& cluster : m_cluster) {
            auto score = m_geneNetwork.findSubgraph(cluster.root, args.depth, args.flexSize, nodes);
            clusterScorePermAbs.emplace_back(std::fabs(score));
            if (!m_permCount) {
                cluster.nodes = nodes;
                cluster.score = score;
            }
        }
    }
    ++m_permCount;
    std::clog << "Perm " << m_permCount << ' ' << perm.get() << std::endl;
    return clusterScorePermAbs;
}

double Experiment::scoreNodeList(const GeneNetwork::NodeList& genes,
                                 const std::vector<int>& label, int nLabels) {
    if (args.sumScore || genes.size() < 2) {  // compute sum(score)
        return m_geneNetwork.subgraphScore(genes);
    }
    else {  // compute score(sum)
        std::vector<double> sum(nSamples());
        for (auto gene : genes) {
            if (args.sumSigned && m_gene[gene].scorePerm < 0)
                sum -= m_gene[gene].expression;
            else
                sum += m_gene[gene].expression;
        }
        return computeScore(sum, label, nLabels);
    }
}

void Experiment::setShrinkageFactors() {
    FastDataSet logVar;
    for (auto& data : m_gene)
        if (data.sd > 0)  // constant probes are likely bad so we ignore them
            logVar.insert(2 * std::log(data.sd));
    EmpiricalBayes eb(logVar.mean(), logVar.var(), m_gene.size(),
                      nSamples() - m_mainPhenotype.nLabels());
    std::cout << "Empirical Bayes: df = " << eb.df()
              << " var = " << eb.var() << std::endl;
    for (auto& data : m_gene)
        data.shrinkageFactor = eb.shrinkageFactor(data.sd * data.sd);
}

size_t Experiment::findPhenotype(const std::string& name) const {
    for (size_t i = 0; i < m_phenotype.size(); ++i)
        if (m_phenotype[i].name() == name)
            return i;
    throw Exception("Unknown phenotype " + name);
}

// Read phenotype file.
// Each line describes a phenotype, in the format:
//   name label_1 label_2 ... label_m

void Experiment::readPhenotypes(const std::string& filename) {
    std::ifstream is(filename.c_str());
    if (!is)
        throw Exception("Could not open " + filename);
    std::string line;
    int lineno = 0;
    while (getline(is, line)) {
        ++lineno;
        std::istringstream ss(line);
        Phenotype phenotype;
        try {
            phenotype.read(ss);
        }
        catch (Exception& e) {
            e << " at " << filename << ':' << lineno;
            throw e;
        }
        m_phenotype.emplace_back(phenotype);
    }

    if (m_phenotype.empty())
        throw Exception("Need at least one phenotype in " + filename);

    // Check names are unique and number of samples are equal.
    std::unordered_set<std::string> names;
    for (auto& p : m_phenotype) {
        auto& name = p.name();
        if (names.find(name) != names.end())
            throw Exception("Phenotype " + name + ": appears twice in " + filename);
        names.insert(name);

        if (p.nSamples() != nSamples())
            throw Exception("Phenotype ") << name << ": has " << p.nSamples()
                                          << " samples expected " << nSamples();
    }

    std::cout << "Read " << m_phenotype.size() << " phenotypes: ";
    for (auto& p : m_phenotype)
        std::cout << p.name() << ' ';
    std::cout << "and " << nSamples() << " samples from " << filename << std::endl;

    size_t ix = args.test.empty() ? 0 : findPhenotype(args.test);
    m_mainPhenotype = m_phenotype[ix];
    if (args.filter.size())
        m_mainPhenotype.filter(args.filter);
}

// Read microarray annotation file.
// Each line describes a probe, in the format:
//     probe_id gene_id
// Gene id is typically the NCBI Gene database id (formerly LocusLink / EntrezGene).
// Gene id "NA" means the probe does not map to a gene; such probes are ignored.

void Experiment::readProbes(const std::string& filename) {
    std::ifstream is(filename.c_str());
    if (!is)
        throw Exception("Could not open " + filename);
    std::string line;
    while (getline(is, line)) {
        std::istringstream ss(line);
        std::string probe, id = "NA";
        ss >> probe >> id;
        if (!ss) {
            throw Exception("Bad probe data " + line);
        }
        else if (id != "NA") {
            auto it = m_probe2geneID.find(probe);
            if (it == m_probe2geneID.end()) {  // new probe
                m_probe2geneID[probe] = id;
                if (m_geneID2index.find(id) == m_geneID2index.end()) {  // new gene
                    m_geneID2index[id] = m_gene.size();
                    GeneData gene;
                    gene.id = id;
                    gene.expression.resize(nSamples());
                    m_gene.emplace_back(gene);
                }
            }
            else if (it->second != id) {
                throw Exception("Probe " + probe + " maps to multiple genes");
            }
        }
    }
    std::cout << "Read " << m_probe2geneID.size() << " probes and "
              << m_gene.size() << " genes from " << filename << std::endl;
}

void Experiment::readGeneNames(const std::string& filename) {
    std::ifstream is(filename.c_str());
    std::string line;
    size_t n = 0;
    while (getline(is, line)) {
        std::istringstream ss(line);
        std::string id, name;
        if (ss >> id >> name) {
            auto it = m_geneID2index.find(id);
            if (it != m_geneID2index.end()) {
                m_gene[it->second].name = name;
                ++n;
            }
        }
    }
    std::cout << "Read " << n << " gene names from " << filename << std::endl;
}

// Read probe expression values from file.
// Each line describes a probe, in the format:
//     probe_id expr_1 expr_2 .. expr_m (one expression value per sample)
// Probes that do not map to a gene (as per the microarray annotation file) are ignored.
// Genes with no probes get expression 0.
// Genes with multiple probes (or multiple copies of the same probe) get probe average.
// Missing / NA expression values are not allowed (for now).

void Experiment::readExpression(const std::string& filename) {
    std::ifstream is(filename.c_str());
    if (!is)
        throw Exception("Could not open " + filename);
    std::string line;
    size_t nProbesRead = 0, nGenesRead = 0;
    while (getline(is, line)) {
        std::istringstream ss(line);
        std::string probe;
        ss >> probe;
        auto it = m_probe2geneID.find(probe);
        if (ss && it != m_probe2geneID.end()) {  // known probe and gene
            std::vector<double> expression;
            for (size_t i = 0; i < nSamples(); i++) {
                double val;
                if (ss >> val)
                    expression.emplace_back(val);
                else
                    break;
            }
            if (!ss)
                throw Exception("Bad probe expression " + line);
            auto& data = m_gene[m_geneID2index[it->second]];
            data.expression += expression;
            ++data.nProbes;
            ++nProbesRead;
        }
    }
    for (auto& data : m_gene) {
        if (data.nProbes) {
            data.expression /= data.nProbes;
            ++nGenesRead;
        }
    }

    std::cout << "Read " << nProbesRead << " expressions for " << nGenesRead
              << " genes from " << filename << std::endl;
}

// Output functions

void Experiment::writeHTML(const std::string& htmlFilename,
                           const std::string& frameFilename,
                           const std::string& startingFrame) const {
    std::cout << "Writing HTML to " << htmlFilename << std::endl;
    std::ofstream os(htmlFilename.c_str());
    os << "<html>" << '\n';
    os << "<title>" << "GXNA " << args.name << ' ' << args.version << "</title>" << '\n';
    os << "<frameset cols=\"35%,65%\">" << '\n';
    os << "<frame src=\"" << frameFilename << "\">" << '\n';
    os << "<frame src=\"" << startingFrame << "\" name=\"frame2\">" << '\n';
    os << "</frameset>" << '\n';
    os << "</html>" << '\n';
}

static void beginHTMLFrame(std::ostream& os) {
    os << "<html>" << '\n';
    os << "<table border=\"1\">" << '\n';
    // Column names
    os << "<tr>" << '\n';
    os << "<th>rank</th>" << '\n';
    os << "<th>root</th>" << '\n';
    os << "<th>rootid</th>" << '\n';
    os << "<th>size</th>" << '\n';
    os << "<th>score</th>" << '\n';
    os << "<th>rawp</th>" << '\n';
    os << "<th>adjp</th>" << '\n';
    os << "</tr>" << '\n';
    os.precision(3);
}

static void endHTMLFrame(std::ostream& os) {
    os << "</table>" << '\n';
    os << "</html>" << '\n';
}

template<typename T>
void addCell(std::ostream& os, const T& val) {
    os << "<td>" << val << "</td>" << '\n';
}

// Add table row for one root/cluster.
static void addRow(std::ostream& os, const std::string& url, size_t n,
                   const std::string& root, const std::string& rootid,
                   size_t size, double score, double rawp, double adjp) {
    os << "<tr>" << '\n';
    if (url.size())
        addCell(os, "<a href=\"" + url + "\" target=\"frame2\">" + std::to_string(n) + "</a>");
    else  // no href
        addCell(os, n);
    addCell(os, root);
    addCell(os, rootid);
    addCell(os, size);
    addCell(os, score);
    addCell(os, rawp);
    addCell(os, adjp);
    os << "</tr>" << '\n';
}

void Experiment::printResults(const MultipleTest<Experiment>& mt, const std::string& path) {
    std::string htmlFilename = "index.html";
    std::string frameFilename = "frame1.html";
    std::string startingFrame;

    std::ofstream osArgs((path + "/" + "parameters.txt").c_str());
    std::ofstream osResults((path + "/" + "results.txt").c_str());
    std::ofstream osFrame((path + "/" + frameFilename).c_str());

    for (size_t i = 0; i < m_gene.size(); ++i)  // set text to draw nodes
        m_geneNetwork.setText(i, m_gene[i].text());

    args.print(osArgs);
    beginHTMLFrame(osFrame);

    std::vector<bool> printed(m_gene.size());
    osResults << std::fixed << std::setprecision(3);
    int nDetailed = 0;
    for (size_t i = 0; i < m_cluster.size(); ++i) {
        auto& cluster = m_cluster[mt.getRank(i)];
        auto& rootData = m_gene[cluster.root];
        auto size = cluster.nodes.size();
        double score = cluster.score;
        double rawP = mt.getRawP(i);
        double adjP = mt.getAdjP(i);

        // Write line to results file.
        osResults << std::setw(5) << i << ' ';
        rootData.printNameId(osResults);
        osResults << std::setw(3) << size << ' '
                  << std::setw(6) << score << ' '
                  << rawP << ' '
                  << adjP << ' ';
        for (auto& v : cluster.nodes)
            osResults << std::setw(6) << m_gene[v].id << ' ';
        osResults << '\n';

        size_t nPrinted = 0;  // check overlap between cluster and previous ones
        for (auto& v : cluster.nodes)
            if (printed[v])
                ++nPrinted;
        if (nPrinted <= args.maxOverlap * size) {  // low overlap, OK to print
            for (auto& v : cluster.nodes)
                printed[v] = true;  // mark as printed
            std::string url;
            if (nDetailed < args.nDetailed) {  // write details
                std::string prefix = "graph_" + std::to_string(nDetailed) + ".";
                std::string filenameTXT = path + "/" + prefix + "txt";
                std::string filenameDOT = path + "/" + prefix + "dot";
                std::string filenameSVG = path + "/" + prefix + "svg";
                bool draw = args.draw && size > 1;  // do not draw singletons
                url = prefix + (draw ? "svg" : "txt");
                if (!startingFrame.size())
                    startingFrame = url;

                m_geneNetwork.write(cluster.nodes, filenameDOT, draw ? filenameSVG : "");
                std::ofstream osTXT(filenameTXT.c_str());
                for (auto& v : cluster.nodes)
                    osTXT << m_gene[v] << '\n';
            }
            if (nDetailed < args.nRows)
                addRow(osFrame, url, i, rootData.name, rootData.id, size, score, rawP, adjP);
            ++nDetailed;
        }
    }
    endHTMLFrame(osFrame);
    writeHTML(path + "/" + htmlFilename, frameFilename, startingFrame);
}

}  // namespace gxna
