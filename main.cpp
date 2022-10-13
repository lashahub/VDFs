#include <iostream>
#include "NTL/ZZ.h"

using namespace NTL;

bool isPrime(const ZZ &n, int k) {

    if (n <= 1) return false;
    if (n <= 3) return true;

    ZZ d = n - 1;
    int s = 0;
    while ((d ^ 1) == 1) {
        d >>= 1;
        s++;
    }

    for (int i = 0; i < k; i++) {
        ZZ a = RandomBits_ZZ(1024) % (n - 3) + 2;
        ZZ x = PowerMod(a, d, n);
        if (x == 1 || x == n - 1) continue;
        for (int j = 0; j < s - 1; j++) {
            x = PowerMod(x, 2, n);
            if (x == n - 1) goto end;
        }
        return false;
        end:;
    }
    return true;

}

int main() {
    int c = 0;
    for (int i = 2; i < 100000; i++) {
        ZZ a(i);
        if (isPrime(a, 100) != isPrime(a, 100)) c++;
    }
    std::cout << c << std::endl;
    return 0;
}
