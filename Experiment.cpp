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

void GeneData::print(std::ostream& os) const {
    os << std::fixed << std::setprecision(3)
       << std::setw(10) << std::left << name << ' '
       << std::setw(6) << std::right << id << ' '
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

inline std::string make_path(const std::string& dir, const std::string& name) {
    return dir + "/" + name;
}

Experiment::Experiment(Args& args_)
    : args(args_)
{
    m_phenotype.read(make_path(args.inputDir, args.phenotypeFile));
    if (args.phenotypes.size())
        m_phenotype.filter(args.phenotypes);
    if (args.useInvariantPerms) { // read secondary phenotype, if needed
        m_phenotypeInvariant.read(make_path(args.inputDir, args.typeFile));
        if (m_phenotypeInvariant.nSamples() != m_phenotype.nSamples())
            throw Exception("Bad number of samples in " + args.typeFile);
    }

    readProbes(make_path(args.refDir, args.probeFile));
    readGeneNames(make_path(args.refDir, args.geneFile));
    readExpression(make_path(args.inputDir, args.expressionFile));
    m_geneNetwork.setSize(m_gene.size()); // so genes with no interactions are still included in the network
    m_geneNetwork.readInteractions(make_path(args.refDir, args.interactionFile), m_geneID2index);

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
    if (args.useInvariantPerms) {
        InvariantPermutationGenerator ppg(m_phenotype.nSamples(), args.nPerms, m_phenotypeInvariant.type());
        ppg.setVerbose(true);
        mt.test(ppg, args.maxTscaled);
    }
    else {
        PermutationGenerator ppg(m_phenotype.nSamples(), args.nPerms);
        ppg.setVerbose(true);
        mt.test(ppg, args.maxTscaled);
    }
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

std::string Experiment::graphFilename(int i, const std::string& type) const {
    return args.name + "_" + args.version + "_" + std::to_string(i) + "." + type; 
}

void Experiment::printResults(const MultipleTest<Experiment>& mt) {
    for (size_t i = 0; i < m_gene.size(); ++i) // need to compute scores first
        m_geneNetwork.setLabel(i, m_gene[i].label());

    auto prefix = args.name + "_" + args.version;
    std::ofstream osArgs(make_path(args.outputDir, prefix + ".arg").c_str());
    std::ofstream osResults(make_path(args.outputDir, prefix + ".res").c_str());
    std::ofstream osHTML(make_path(args.outputDir, prefix + ".html").c_str());
    auto frameFilename = prefix + "_frame1.html";
    std::ofstream osFrame(make_path(args.outputDir, frameFilename).c_str());

    args.print(osArgs);
    beginHTML(osHTML, osFrame, frameFilename); // prepare html files
    std::vector<bool> printed(m_gene.size());
    osResults.precision(4);
    int nDOT = 0;
    for (int i = 0; i < m_testData.size(); ++i) {
        auto& testData = m_testData[mt.getRank(i)];
        int root = testData.root;
        osResults << i << ' ' << m_gene[root] << ' ';
        for (auto& v : testData.cluster)
            osResults << m_gene[v].id << ' ';
        double rawP = mt.getRawP(i), adjP = mt.getAdjP(i);
        osResults << testData.score << ' ' << rawP << ' ' << adjP << '\n';

        // now check if we print to html and make graphs
        int n1 = 0, n2 = testData.cluster.size();
        for (auto& v : testData.cluster)
            if (printed[v])
                ++n1;
        if (n1 <= args.maxOverlap * n2) { // ok to print
            for (auto& v : testData.cluster)
                printed[v] = true; // mark as printed
            if (nDOT < args.graphCount) {
                std::string filenameDOT = make_path(args.outputDir, graphFilename(nDOT, "dot"));
                std::string filenameSVG = make_path(args.outputDir, graphFilename(nDOT, "svg"));
                m_geneNetwork.write(testData.cluster, filenameDOT, args.draw ? filenameSVG : "");
                std::ofstream osTXT(make_path(args.outputDir, graphFilename(nDOT, "txt")).c_str());
                for (auto& v : testData.cluster)
                    osTXT << m_gene[v] << '\n';
            }
            if (nDOT < args.maxRows)
                addRow(osFrame, i, nDOT, m_gene[root].name, m_gene[root].id, n2, testData.score, rawP, adjP);
            ++nDOT;
        }
    }
    endHTML(osHTML);
}

void Experiment::beginHTML(std::ostream& osHTML, std::ostream& osFrame, const std::string& frameFilename) const {
    std::string startingFrame = graphFilename(0, args.draw ? "svg" : "txt");

    // write main html file
    osHTML << "<html>" << '\n';
    osHTML << "<title>" << "GXNA " << args.name << ' ' << args.version << "</title>" << '\n';
    osHTML << "<frameset cols=\"35%,65%\">" << '\n';
    osHTML << "<frame src=\"" << frameFilename << "\">" << '\n';
    osHTML << "<frame src=\"" << startingFrame << "\" name=\"frame2\">" << '\n';
    osHTML << "</frameset>" << '\n';
    osHTML << "</html>" << '\n';

    // start frame html file
    osFrame << "<html>" << '\n';
    osFrame << "<table border=\"1\">" << '\n';
    osFrame << "<tr>" << '\n';
    osFrame << "<th>rank</th>" << '\n';
    osFrame << "<th>root</th>" << '\n';
    osFrame << "<th>rootid</th>" << '\n';
    osFrame << "<th>size</th>" << '\n';
    osFrame << "<th>score</th>" << '\n';
    osFrame << "<th>rawp</th>" << '\n';
    osFrame << "<th>adjp</th>" << '\n';
    osFrame << "</tr>" << '\n';
    osFrame.precision(3);
}

void Experiment::endHTML(std::ostream& os) const {
    os << "</table>" << '\n';
    os << "</html>" << '\n';
}
void Experiment::addRow(std::ostream& os, int n, int nDOT, const std::string& root, const std::string& rootid, int size,
                        double score, double rawp, double adjp) const {                 
    os << "<tr>" << '\n';
    if (nDOT < args.graphCount) {
        std::string url = graphFilename(nDOT, args.draw ? "svg" : "txt");
        os << "<td><a href=\"" << url << "\" target=\"frame2\">" << n << "</a></td>" << '\n';
    }
    else
        os << "<td>" << n << "</td>" << '\n';
    os << "<td>" << root << "</td>" << '\n';
    os << "<td>" << rootid << "</td>" << '\n';
    os << "<td>" << size << "</td>" << '\n';
    os << "<td>" << score << "</td>" << '\n';
    os << "<td>" << rawp << "</td>" << '\n';
    os << "<td>" << adjp << "</td>" << '\n';
    os << "</tr>" << '\n';
}

} // namespace gxna
