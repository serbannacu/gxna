#include "Experiment.h"
#include "Distribution.h"
#include "Exception.h"
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
        if (it == m_label2n.end()) {
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

std::string GeneData::label() const {
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

    for (auto& data : m_gene) {
        FastDataSet d(data.expression);
        data.sd = d.sd();
    }

    if (args.shrink)  // use empirical Bayes shrinkage
        setShrinkageFactor();
}

void Experiment::run() {
    // Select gene roots for multiple testing
    // Potential roots will not be within args.minDistance of each other

    std::vector<bool> ignore(m_gene.size());
    for (size_t root = 0; root < m_geneNetwork.nNodes(); ++root)
        if (m_geneNetwork.degree(root) >= unsigned(args.minDegree)
            && m_gene[root].sd >= args.minSD && !ignore[root]) {
            TestData testData;
            testData.root = root;
            if (args.algoType == AlgoType::Basic)
                testData.cluster = m_geneNetwork.ball(root, args.radius);
            if (args.minDistance > 0) {  // do not use nearby genes as roots
                for (auto v : m_geneNetwork.ball(root, args.minDistance))
                    ignore[v] = true;
            }
            m_testData.emplace_back(testData);
        }

    // Prepare permutations for resampling

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
    MultipleTest<Experiment> mt(*this, m_testData.size(), Eps);
    std::cout << "Testing " << m_testData.size() << " objects using phenotype "
              << m_mainPhenotype << std::endl;
    mt.test(*pg, args.maxTscaled);
    delete pg;

    for (size_t i = 0; i < m_gene.size(); ++i)  // set labels before calling printResults
        m_geneNetwork.setLabel(i, m_gene[i].label());  // label includes gene score
    printResults(mt);
}

static double computeScore(const std::vector<double>& expression,
                           const std::vector<int>& pheno, int nLabels) {
    double score;
    if (nLabels > 2) {  // F statistic
        double f = fstatPheno(&expression[0], pheno, nLabels);
        score = f2z(f, nLabels - 1, pheno.size() - nLabels);
    }
    else {  // t statistic; may want to convert to z-score
        double t = tstatPheno(&expression[0], pheno, 0, 1);  // pheno 0 vs 1
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
    // Recompute gene scores

    auto pheno = perm.apply(m_mainPhenotype.get());
    auto nLabels = m_mainPhenotype.nLabels();

    std::vector<double> geneScorePermAbs;
    for (auto& data : m_gene) {
        double score = computeScore(data.expression, pheno, nLabels);
        if (args.shrink)
            score *= data.shrinkageFactor;
        data.scorePerm = score;
        if (!m_permCount)
            data.score = score;
        geneScorePermAbs.emplace_back(std::fabs(score));
    }
    m_geneNetwork.setScores(geneScorePermAbs, args.scalingExponent);

    // Recompute cluster scores

    std::vector<double> clusterScorePermAbs;
    clusterScorePermAbs.reserve(m_testData.size());
    if (args.algoType == AlgoType::Basic) {
        for (auto& testData : m_testData) {
            auto val = scoreNodeList(testData.cluster, pheno, nLabels);
            clusterScorePermAbs.emplace_back(std::fabs(val));
            if (!m_permCount)
                testData.score = val;
        }
    }
    else {  // AlgoType::GXNA
        GeneNetwork::NodeList cluster;
        for (auto& testData : m_testData) {
            auto val = m_geneNetwork.findSubgraph(testData.root, args.depth, args.flexSize, cluster);
            clusterScorePermAbs.emplace_back(std::fabs(val));
            if (!m_permCount) {
                testData.cluster = cluster;
                testData.score = val;
            }
        }
    }
    ++m_permCount;
    return clusterScorePermAbs;
}

// Scoring function for MultipleTest, predefined case
double Experiment::scoreNodeList(const GeneNetwork::NodeList& genes,
                                 const std::vector<int>& pheno, int nLabels) {
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
        return computeScore(sum, pheno, nLabels);
    }
}

void Experiment::setShrinkageFactor() {
    FastDataSet logVar;
    for (auto& data : m_gene)
        if (data.sd > 0)  // constant probes are likely bad so we ignore them
            logVar.insert(2 * std::log(data.sd));
    EmpiricalBayes eb;
    eb.estimate(logVar.mean(), logVar.var(), m_gene.size(),
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
//   name sample_1 sample_2 ... sample_n

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

// Read probe expression values from file
// Each line describes a probe, in the format:
//     probe_id expr_1 expr_2 .. expr_n (one expression per phenotype)
// Probes that do not map to a gene (as per the microarray annotation file) are ignored
// Genes with no probes get expression 0
// Genes with multiple probes (or multiple copies of the same probe) get probe average
// Missing / NA expression values are not allowed (for now)

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
        if (ss && it != m_probe2geneID.end()) {
            std::vector<double> expression;
            for (size_t i = 0; i < nSamples(); i++) {
                double val;
                if (ss >> val)
                    expression.emplace_back(val);
                else
                    break;
            }
            if (ss) {
                auto& data = m_gene[m_geneID2index[it->second]];
                ++data.nProbes;
                data.expression += expression;
                ++nProbesRead;
            }
            else
                throw Exception("Bad probe expression " + line);
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

static void addRow(std::ostream& os, const std::string& url, size_t n,
                   const std::string& root, const std::string& rootid,
                   size_t size, double score, double rawp, double adjp) {
    os << "<tr>" << '\n';
    if (url.size())
        addCell(os, "<a href=\"" + url + "\" target=\"frame2\">" + std::to_string(n) + "</a>");
    else
        addCell(os, n);
    addCell(os, root);
    addCell(os, rootid);
    addCell(os, size);
    addCell(os, score);
    addCell(os, rawp);
    addCell(os, adjp);
    os << "</tr>" << '\n';
}

void Experiment::printResults(const MultipleTest<Experiment>& mt) {
    std::string path = args.outputDir + "/" + args.name + "/" + args.version;
    std::filesystem::create_directories(path);

    std::string htmlFilename = "index.html";
    std::string frameFilename = "frame1.html";
    std::string startingFrame;

    std::ofstream osArgs((path + "/" + "parameters.txt").c_str());
    std::ofstream osResults((path + "/" + "results.txt").c_str());
    std::ofstream osFrame((path + "/" + frameFilename).c_str());

    args.print(osArgs);
    beginHTMLFrame(osFrame);

    std::vector<bool> printed(m_gene.size());
    osResults << std::fixed << std::setprecision(3);
    int nDetailed = 0;
    for (size_t i = 0; i < m_testData.size(); ++i) {
        auto& testData = m_testData[mt.getRank(i)];
        auto& rootData = m_gene[testData.root];
        auto size = testData.cluster.size();
        double score = testData.score, rawP = mt.getRawP(i), adjP = mt.getAdjP(i);

        osResults << std::setw(5) << i << ' ';
        rootData.printNameId(osResults);
        osResults << std::setw(3) << size << ' '
                  << std::setw(6) << score << ' '
                  << rawP << ' '
                  << adjP << ' ';
        for (auto& v : testData.cluster)
            osResults << std::setw(6) << m_gene[v].id << ' ';
        osResults << '\n';

        // Now check overlap between current and previous clusters
        size_t nPrinted = 0;
        for (auto& v : testData.cluster)
            if (printed[v])
                ++nPrinted;
        if (nPrinted <= args.maxOverlap * size) {  // low overlap, OK to print
            for (auto& v : testData.cluster)
                printed[v] = true;  // mark as printed
            std::string url;
            if (nDetailed < args.nDetailed) {
                std::string prefix = "graph_" + std::to_string(nDetailed) + ".";
                std::string filenameTXT = path + "/" + prefix + "txt";
                std::string filenameDOT = path + "/" + prefix + "dot";
                std::string filenameSVG = path + "/" + prefix + "svg";
                bool draw = args.draw && size > 1;  // do not draw singletons
                url = prefix + (draw ? "svg" : "txt");
                if (!startingFrame.size())
                    startingFrame = url;

                m_geneNetwork.write(testData.cluster, filenameDOT, draw ? filenameSVG : "");
                std::ofstream osTXT(filenameTXT.c_str());
                for (auto& v : testData.cluster)
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
