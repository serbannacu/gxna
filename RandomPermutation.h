// Random permutation generators

#pragma once

#include <iosfwd>
#include <vector>

namespace gxna {

// permutation of the integers {0, 1, 2, ..., n - 1}
class Permutation {
public:
    Permutation(size_t n)
        : m_v(n)
    {
        init();
    }

    const std::vector<int>& get() const { return m_v; }

    int operator[](size_t i) const { return m_v[i]; }

    template<typename T>
    std::vector<T> apply(const std::vector<T>& v) const {
        std::vector<T> w(v.size());
        for (size_t i = 0; i < v.size(); ++i)
            w[i] = v[m_v[i]];
        return w;
    }

    void init() { // to Id permutation
        for (size_t i = 0; i < m_v.size(); ++i)
            m_v[i] = i;
    }
    
    void randomize(size_t start, size_t len); // uniform random scramble obtained via element swaps
    
    void randomize() {
        randomize(0, m_v.size());
    }

private:
    std::vector<int> m_v;
};

class PermutationHistogram {
public:
    PermutationHistogram(size_t n) {
        m_count.resize(n);
        for (auto& v : m_count)
            v.resize(n);
        m_n = n;
        m_n_add = 0;
    }

    void add(const Permutation& p) {
        for (size_t i = 0; i < m_n; ++i)
            ++m_count[i][p[i]];
        ++m_n_add;
    }

    void print(std::ostream& os, int prec = 4) const;
    
private:
    std::vector< std::vector<int> > m_count;
    size_t m_n;
    size_t m_n_add;
};

// PermutationGenerator yields a sequence of uniform random permutations of the integers {0, 1, 2, ..., n - 1}
// first permutation returned is guaranteed to be Id

class PermutationGenerator {
public:
    PermutationGenerator(size_t n, size_t limit)
        : m_perm(n),
          m_count(0),
          m_limit(limit),
          m_percentage(0),
          m_verbose(false)
    {}

    void setVerbose(bool verbose = true) {
        m_verbose = verbose;
    }

    size_t count() const {
        return m_count;
    }
    
    const Permutation& get() const {
        return m_perm;
    }
    
    bool next();

protected:
    Permutation m_perm;

private:
    virtual void update() {
        m_perm.randomize();
    }

    void showProgress();

    size_t m_count; // number of permutations already generated
    size_t m_limit; // total number of permutations to generate
    int m_percentage;
    bool m_verbose;
};

// InvariantPermutationGenerator yields a sequence of uniform invariant random permutations
// first permutation returned is guaranteed to be Id
// Each element in {0, 1, 2, ..., n - 1} has an assigned type
// Invariance means all permutations $p$ preserve type, so $p(i)$ has the same type as $i$ for all $i$
// Thus each permutation generated induces an uniform random permutation within each type

// the ctor expectes a vector of types in sequential order, so {0, 0, 1, 1, 1} is OK but {0, 1, 0, 1, 1} is not

class InvariantPermutationGenerator : public PermutationGenerator {
public:
    InvariantPermutationGenerator(size_t n, size_t limit, const std::vector<int>& types)
        : PermutationGenerator(n, limit)
    {
        assert(n == types.size());
        int type;
        for (auto val : types) {
            if (m_typeCount.empty() || type != val) {
                type = val;
                m_typeCount.emplace_back(0);
            }
            ++m_typeCount.back();
        }
    }

private:
    void update() {
        size_t sum = 0;
        for (auto val : m_typeCount) {
            m_perm.randomize(sum, val);
            sum += val;
        }
    }

    std::vector<int> m_typeCount;
};

} // namespace gxna
