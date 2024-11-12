#pragma once

#include <iosfwd>
#include <vector>

namespace gxna {

// Permutation of the integers {0, 1, 2, ..., n - 1}

class Permutation {
 public:
    explicit Permutation(size_t n)
        : m_v(n) {
        init();
    }

    const std::vector<int>& get() const {
        return m_v;
    }

    size_t size() const {
        return m_v.size();
    }

    int operator[](size_t i) const {
        return m_v[i];
    }

    // Seed random number generator
    static void seed(int);

    // Scramble vector
    template<typename T>
    std::vector<T> apply(const std::vector<T>& v) const {
        std::vector<T> w(v.size());
        for (size_t i = 0; i < v.size(); ++i)
            w[i] = v[m_v[i]];
        return w;
    }

    // Initialize to Id permutation
    void init() {
        for (size_t i = 0; i < m_v.size(); ++i)
            m_v[i] = i;
    }

    // Uniform random scramble obtained via element swaps
    void randomize(size_t start, size_t len);

    void randomize() {
        randomize(0, m_v.size());
    }

 private:
    std::vector<int> m_v;
};

// For a set of permutations, PermutationHistogram tracks the distribution
// of each position. It can be used to check if a sample of permutations
// is approximately uniformly distributed.

class PermutationHistogram {
 public:
    explicit PermutationHistogram(size_t n)
        : m_count(n),
          m_n(n),
          m_n_perm(0) {
        for (auto& v : m_count)
            v.resize(n);
    }

    const auto& getCount(size_t i) {
        return m_count[i];
    }

    // Add permutation to histogram.
    void insert(const Permutation& p);

    // For each position, print the distribution of its value.
    void print(std::ostream& os, int prec = 4) const;

 private:
    std::vector< std::vector<int> > m_count;
    size_t m_n;  // size of each permutation
    size_t m_n_perm;  // number of permutations
};

// PermutationGenerator yields a sequence of permutations
// of the integers {0, 1, 2, ..., n - 1}.
// The first permutation returned is guaranteed to be Id.

class PermutationGenerator {
 public:
    PermutationGenerator(size_t n, size_t limit)
        : m_perm(n),
          m_count(0),
          m_limit(limit),
          m_percentage(0),
          m_verbose(false)
    {}

    virtual ~PermutationGenerator() = default;

    void setVerbose(bool verbose = true) {
        m_verbose = verbose;
    }

    size_t count() const {
        return m_count;
    }

    const Permutation& get() const {
        return m_perm;
    }

    // Generate next permutation, return false if the end is reached.
    bool next();

 protected:
    Permutation m_perm;

 private:
    // Generate next permutation, store in m_perm.
    virtual void update() = 0;

    void showProgress();

    size_t m_count;  // number of permutations already generated
    size_t m_limit;  // total number of permutations to generate
    int m_percentage;  // used to show progress
    bool m_verbose;
};

// Uniform random permutation

class UniformPermutation : public PermutationGenerator {
 public:
    UniformPermutation(size_t n, size_t limit)
        : PermutationGenerator(n, limit)
    {}

 private:
    virtual void update() {
        m_perm.randomize();
    }
};

// InvariantPermutation yields a sequence of uniform INVARIANT random permutations.
// Each element in {0, 1, 2, ..., n - 1} has an assigned label.
// Invariance means all permutations $p$ preserve labels,
// so $p(i)$ has the same label as $i$ for all $i$.
// Thus each invariant permutation induces an uniform random permutation within each type.

// The ctor expectes a vector of labels in sequential order,
// so {0, 0, 1, 1, 1} is OK but {0, 1, 0, 1, 1} is not.

class InvariantPermutation : public PermutationGenerator {
 public:
    InvariantPermutation(size_t n, size_t limit, const std::vector<int>& label);

 private:
    virtual void update();

    std::vector<int> m_labelCount;
};

}  // namespace gxna
