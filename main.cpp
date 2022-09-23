#include <iostream>
#include "NTL/ZZ.h"

using namespace NTL;

int main() {
    ZZ a, b, c;
    std::cin >> a >> b;
    c = (a + 1) * (b + 1);
    std::cout << c << std::endl;
    return 0;
}
