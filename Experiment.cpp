#include "Experiment.h"
#include "Distribution.h"
#include "Exception.h"
#include "MultipleTest.h"
#include "RandomPermutation.h"
#include "Statistics.h"
#include "VectorUtil.h"

#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <filesystem>

namespace gxna {

void PhenotypeVector::read(const std::string& filename) {
    std::ifstream is(filename.c_str());
    if (!is)
        throw Exception("Could not open " + filename);
    std::string name;
    while (is >> name) {
        auto it = m_name2type.find(name);
        int i;
        if (it == m_name2type.end()) {
            m_name2type[name] = i = m_type2name.size();
            m_type2name.emplace_back(name);
        }
        else
            i = it->second;
        m_type.emplace_back(i);
    }
    if (nTypes() < 2)
        throw Exception("Need at least two phenotypes in " + filename);

    std::cerr << "Read " << nSamples() << " samples and " << nTypes() << " types from " << filename << ": ";
    for (auto& name : m_type2name)
        std::cerr << name << ' ';
    std::cerr << '\n';
}

void PhenotypeVector::filter(const std::vector<std::string>& v) {
    for (auto& name : v) {
        auto it = m_name2type.find(name);
        if (it == m_name2type.end())
            throw Exception("Unknown phenotype " + name);
        // should use it->second
    }
    throw Exception("PhenotypeVector::filter not implemented");
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
       << std::setw(7) << score
        ;
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
    : args(args_)
{
    m_phenotype.read(args.inputDir + "/" + args.phenotypeFile);
    if (args.phenotypes.size())
        m_phenotype.filter(args.phenotypes);
    if (args.invariantPerms) { // read secondary phenotype, if needed
        m_phenotypeInvariant.read(args.inputDir + "/" + args.typeFile);
        if (m_phenotypeInvariant.nSamples() != m_phenotype.nSamples())
            throw Exception("Bad number of samples in " + args.typeFile);
    }

    readProbes(args.refDir + "/" + args.probeFile);
    readGeneNames(args.refDir + "/" + args.geneFile);
    readExpression(args.inputDir + "/" + args.expressionFile);
    m_geneNetwork.setSize(m_gene.size()); // so genes with no interactions are still included in the network
    m_geneNetwork.readInteractions(args.refDir + "/" + args.interactionFile, m_geneID2index);

    for (auto& data : m_gene) {
        FastDataSet d(data.expression);
        data.sd = d.sd();
    }
    if (args.shrink) // use empirical Bayes shrinkage
        setShrinkageFactor();
}

void Experiment::run() {
    // select gene roots for multiple testing
    // for AlgoType::GXNA potential roots will not be within args.radius of each other
    std::vector<bool> ignore(m_gene.size());
    for (int root = 0; root < m_geneNetwork.nNodes(); ++root)
        if (m_geneNetwork.degree(root) >= args.minDegree && m_gene[root].sd >= args.minSD && !ignore[root]) {
            auto ball = m_geneNetwork.ball(root, args.radius);
            TestData testData;
            testData.root = root;
            if (args.algoType == AlgoType::Basic) {
                testData.cluster = ball;
            }
            else { // AlgoType::GXNA
                for (auto v : ball)
                    ignore[v] = true; // do not use nearby genes as roots
            }
            m_testData.emplace_back(testData);
        }
    std::cerr << "Testing " << m_testData.size() << " objects\n";
    
    MultipleTest<Experiment> mt(*this, m_testData.size());
    PermutationGenerator *pg;
    if (args.invariantPerms)
        pg = new InvariantPermutationGenerator(m_phenotype.nSamples(), args.nPerms, m_phenotypeInvariant.type());
    else
        pg = new UniformPermutationGenerator(m_phenotype.nSamples(), args.nPerms);
    pg->setVerbose(true);
    mt.test(*pg, args.maxTscaled);
    delete pg;

    for (size_t i = 0; i < m_gene.size(); ++i) // need to set labels before calling printResults
        m_geneNetwork.setLabel(i, m_gene[i].label()); // label includes gene score

    printResults(mt);
}

static double computeScore(const std::vector<double>& expression, const std::vector<int>& phenotype, int nPhenotypes) {
    if (nPhenotypes > 2) { // F statistic
        double f = fstat(&expression[0], phenotype, nPhenotypes);
        return zfCDF(f, nPhenotypes - 1, phenotype.size() - nPhenotypes);
    }
    else { // T statistic; may want to convert to z-score
        double t = tstatPheno(&expression[0], phenotype, 0, 1); // pheno 0 vs 1
        return t;
    }
}

std::vector<double> Experiment::operator()(const Permutation& perm) {
    // Recompute gene scores
    auto phenotype = perm.apply(m_phenotype.type());
    std::vector<double> geneScorePermAbs;
    for (auto& data : m_gene) {
        double score = computeScore(data.expression, phenotype, m_phenotype.nTypes());
        if (args.shrink)
            score *= data.shrinkageFactor;
        data.scorePerm = score;
        if (!m_permCount)
            data.score = score;
        geneScorePermAbs.emplace_back(fabs(score));
    }
    m_geneNetwork.setScores(geneScorePermAbs, args.scalingExponent);
  
    // Recompute cluster scores
    std::vector<double> clusterScorePermAbs;
    clusterScorePermAbs.reserve(m_testData.size());
    if (args.algoType == AlgoType::Basic) {
        for (auto& testData : m_testData) {
            auto val = scoreNodeList(testData.cluster, phenotype);
            clusterScorePermAbs.emplace_back(fabs(val));
            if (!m_permCount)
                testData.score = val;
        }
    }
    else { // AlgoType::GXNA
        std::vector<int> cluster;
        for (auto& testData : m_testData) {
            auto val = m_geneNetwork.findSubgraph(testData.root, args.depth, args.flexSize, cluster);
            clusterScorePermAbs.emplace_back(fabs(val));
            if (!m_permCount) {
                testData.cluster = cluster;
                testData.score = val;
            }
        }
    }
    ++m_permCount;
    return clusterScorePermAbs;
}

void Experiment::setShrinkageFactor() {
    FastDataSet logVar;
    for (auto& data : m_gene)
        if (data.sd > 0) // constant probes are presumably bad so we ignore them
            logVar.insert(2 * log(data.sd));
    EmpiricalBayes eb;
    eb.estimate(logVar.mean(), logVar.var(), m_gene.size(), m_phenotype.nSamples() - m_phenotype.nTypes());
    std::cerr << "Shrinkage parameters: df = " << eb.df() << " var = " << eb.var() << '\n';
    for (auto& data : m_gene)
        data.shrinkageFactor = eb.shrinkageFactor(data.sd * data.sd);
}

// scoring function for MultipleTest, predefined case
double Experiment::scoreNodeList(const GeneNetwork::NodeList& genes, const std::vector<int>& phenotype) {
    if (args.sumScore || genes.size() < 2) { // compute sum(score)
        return m_geneNetwork.subgraphScore(genes);
    }
    else { // compute score(sum)
        std::vector<double> sum(m_phenotype.nSamples());
        for (int gene : genes) {
            if (args.sumSigned && m_gene[gene].scorePerm < 0)
                sum -= m_gene[gene].expression;
            else
                sum += m_gene[gene].expression;
        }
        return computeScore(sum, phenotype, m_phenotype.nTypes());
    }
}

// Read microarray annotation file
// Each line describes a probe, in the format: probe_id gene_id
// Gene id is typically the NCBI Gene database id (formerly LocusLink / EntrezGene)
// Gene id "NA" means the probe does not map to a gene; such probes are ignored

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
            std::cerr << "readProbes bad line " << line << '\n';
        }
        else if (id != "NA") {
            auto it = m_probe2geneID.find(probe);
            if (it == m_probe2geneID.end()) { // new probe
                m_probe2geneID[probe] = id;
                if (m_geneID2index.find(id) == m_geneID2index.end()) { // new gene
                    m_geneID2index[id] = m_gene.size();
                    GeneData gene;
                    gene.id = id;
                    gene.expression.resize(m_phenotype.nSamples());
                    m_gene.emplace_back(gene);
                }
            }
            else if (it->second != id) { 
                std::cerr << "readProbes probe " << probe << " maps to multiple genes "
                          << it->second << ' ' << id << '\n';
            }
        }
    }
    std::cerr << "Read " << m_probe2geneID.size() << " probes and "
              << m_gene.size() << " genes from " << filename << '\n';
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
    std::cerr << "Read " << n << " gene names from " << filename << "\n";
}

// Read probe expression values from file
// Each line describes a probe, in the format: probe_id expr_1 expr_2 .. expr_n (one expression per phenotype)
// Probes that do not map to a gene (as per the microarray annotation file) are ignored
// Genes with no probes get expression 0
// Genes with multiple probes (or multiple copies of the same probe) get probe average (for now)
// Missing / NA expression values are not allowed (for now)

void Experiment::readExpression(const std::string& filename) {
    std::ifstream is(filename.c_str());
    if (!is)
        throw Exception("Could not open " + filename);
    std::string line;
    int nProbesRead = 0, nGenesRead = 0;
    while (getline(is, line)) {
        std::istringstream ss(line);
        std::string probe;
        ss >> probe;
        auto it = m_probe2geneID.find(probe);
        if (ss && it != m_probe2geneID.end()) {
            std::vector<double> expression;
            for (int i = 0; i < m_phenotype.nSamples(); i++) {
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
                std::cerr << "readExpression bad line " << line << '\n';
        }
    }
    for (auto& data : m_gene) {
        if (data.nProbes) {
            data.expression /= data.nProbes;
            ++nGenesRead;
        }
    }

    std::cerr << "Read " << nProbesRead << " expressions for " << nGenesRead << " genes from " << filename << '\n';
}

// output functions

void Experiment::writeHTML(const std::string& htmlFilename, const std::string& frameFilename, const std::string& startingFrame) const {
    std::cerr << "Writing output to " << htmlFilename << '\n';
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

static void addRow(std::ostream& os, const std::string& url, int n, const std::string& root, const std::string& rootid, int size,
                   double score, double rawp, double adjp) {                 
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
    for (int i = 0; i < m_testData.size(); ++i) {
        auto& testData = m_testData[mt.getRank(i)];
        auto& rootData = m_gene[testData.root];
        int size = testData.cluster.size();
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

        // now check overlap between current and previous clusters
        int nPrinted = 0;
        for (auto& v : testData.cluster)
            if (printed[v])
                ++nPrinted;
        if (nPrinted <= args.maxOverlap * size) { // low overlap, OK to print
            for (auto& v : testData.cluster)
                printed[v] = true; // mark as printed
            std::string url;
            if (nDetailed < args.nDetailed) {
                std::string prefix = "graph_" + std::to_string(nDetailed) + ".";
                std::string filenameTXT = path + "/" + prefix + "txt";
                std::string filenameDOT = path + "/" + prefix + "dot";
                std::string filenameSVG = path + "/" + prefix + "svg";
                bool draw = args.draw && size > 1; // do not draw singletons
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

} // namespace gxna
