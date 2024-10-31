#include "Distribution.h"

#include <iostream>

using namespace gxna;

void test_normCDF() {
    for (double x = 0.1; x < 2; x += 0.1)
        std::cout << x << ' ' << normCDF(x) << '\n';
}

int main() {
    test_normCDF();
    return 0;
}
