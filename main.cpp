#include <iostream>
#include <cmath>
#include <chrono>
#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "bicycl.hpp"

using namespace NTL;

class RSA_Group {

public:
    int k;
    ZZ p, q, N, fi;


public:

    RSA_Group(const ZZ &p, const ZZ &q) {
        this->p = p;
        this->q = q;
        N = p * q;
        fi = (p - 1) * (q - 1);
    }

    ZZ mul(const ZZ &a, const ZZ &b) const {
        return MulMod(a, b, N);
    }

    ZZ trapdoor(const ZZ &g, uint64_t t) const {
        ZZ two = ZZ(2);
        ZZ zt = ZZ(t);
        ZZ e = PowerMod(two, zt, fi);
        return PowerMod(g, e, N);
    }

    ZZ sequential(const ZZ &g, uint64_t t) const {
        ZZ y = g;
        ZZ two = ZZ(2);
        for (uint64_t i = 0; i < t; i++) {
            y = PowerMod(y, two, N);
        }
        return y;
    }

    ZZ genproof(const ZZ &g, const ZZ &t, const ZZ &l) const {
        ZZ pow = (2 ^ t) / l;
        return PowerMod(g, pow, N);
    }

};

class Class_Group {

public:

    ZZ generateD() {
        ZZ d = -RandomPrime_ZZ(1024, 50);
        while (d % 4 != 1) {
            d = -RandomPrime_ZZ(1024, 50);
        }
        return d;
    }

    ZZ MinkowskiBound() {
        ZZ delta = SqrRoot(-generateD());
        ZZ bound = 2 * delta / M_PI;
        return bound;
    }

    Class_Group() {

    }

};

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
    std::cout<<"HI"<<std::endl;
    int k = 1024;
    std::cout<<"HI"<<std::endl;
    RSA_Group group(RandomPrime_ZZ(k, 50), RandomPrime_ZZ(k, 50));
    std::cout<<"HI"<<std::endl;
    ZZ g = RandomBits_ZZ(k);
    uint64_t time = 200;

    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << group.trapdoor(g, time) << std::endl << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << ms.count() << std::endl;

    t1 = std::chrono::high_resolution_clock::now();

    ZZ s = group.sequential(g, time);

    std::cout << s << std::endl << std::endl;

    t2 = std::chrono::high_resolution_clock::now();
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << ms.count() << std::endl;

    t1 = std::chrono::high_resolution_clock::now();

    ZZ l = RandomPrime_ZZ(2 * k, 50);

    ZZ pi = group.genproof(g, ZZ(time), l);
    ZZ r = PowerMod(ZZ(2), ZZ(time), l);

    ZZ first = PowerMod(pi, l, group.N);
    ZZ second = PowerMod(g, r, group.N);
    ZZ mul = MulMod(first, second, group.N);

    bool eq = (mul == s);

    std::cout << (eq ? "Verified" : "Failed") << std::endl << std::endl;

    t2 = std::chrono::high_resolution_clock::now();
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << ms.count() << std::endl;

    return 0;
}
