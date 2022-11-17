#include <iostream>
#include <cmath>
#include <chrono>
#include "NTL/ZZ.h"
#include "NTL/ZZ_p.h"
#include "bicycl.hpp"
#include "openssl/evp.h"
#include "openssl/sha.h"

#define NUM_TRIALS 50
#define SECURITY_PARAM 128
#define NUM_LENGTH 1536


using namespace NTL;
using namespace BICYCL;

ZZ to_int(const std::string &x) {
    ZZ res;
    for (char c: x) {
        res *= 2;
        res += (c == '1') ? 1 : 0;
    }
    return res;
}

std::string to_bin(ZZ x) {
    if (x == 0) return "0";
    std::string res;
    while (x > 0) {
        res += (x % 2 == 1) ? '1' : '0';
        x /= 2;
    }
    std::reverse(res.begin(), res.end());
    return res;
}

class RSA_Group {

public:
    ZZ p, q, N, fi;

public:

    RSA_Group() {
        RandomPrime(this->p, NUM_LENGTH, NUM_TRIALS);
        RandomPrime(this->q, NUM_LENGTH, NUM_TRIALS);
        N = p * q;
        ZZ_p::init(N);
        fi = (p - 1) * (q - 1);
    }

    std::string to_str(uint8_t *start, int length) {
        std::string res;
        for (int i = 0; i < length; i++) {
            uint8_t curr = start[i];
            std::string s;
            for (int j = 0; j < 8; j++) {
                s += (curr % 2 == 1) ? '1' : '0';
                curr /= 2;
            }
            std::reverse(s.begin(), s.end());
            res += s;
        }
        return res;
    }

    ZZ hash(std::string x) {

        //std::string bin = to_bin(x);
        std::string bin = "0110011001";
        const char *cstr = bin.c_str();
        const unsigned char *input = reinterpret_cast<unsigned char *>(const_cast<char *>(cstr));
        unsigned char output[32 + 1];
        SHA3_256(output, input, bin.length());
        output[32] = '\0';
        return to_int(std::string(reinterpret_cast<char *>(output)));


        /*std::string input = "residue" + x;
        uint32_t digest_length = SHA256_DIGEST_LENGTH;
        const EVP_MD *algorithm = EVP_sha3_256();
        auto *digest = static_cast<uint8_t *>(OPENSSL_malloc(digest_length));
        EVP_MD_CTX *context = EVP_MD_CTX_new();
        EVP_DigestInit_ex(context, algorithm, nullptr);
        EVP_DigestUpdate(context, input.c_str(), input.size());
        EVP_DigestFinal_ex(context, digest, &digest_length);
        EVP_MD_CTX_destroy(context);
        std::string output = to_str(digest, digest_length);
        OPENSSL_free(digest);
        return to_int(output) % N;*/
    }

    ZZ trapdoor(const ZZ &g, uint64_t t) const {
        ZZ two = ZZ(2);
        ZZ zt = ZZ(t);
        ZZ e = PowerMod(two, zt, fi);
        return PowerMod(g, e, N);
    }

    ZZ eval(const ZZ &g, uint64_t t) const {
        ZZ y = g;
        ZZ two = ZZ(2);
        for (uint64_t i = 0; i < t; i++) {
            y = PowerMod(y, two, N);
        }
        return y;
    }

    ZZ trapdoorproof(const ZZ &g, const ZZ &t, const ZZ &l) {

    }

    ZZ genproof(const ZZ &g, const ZZ &t, const ZZ &l) const {
        ZZ pow = (2 ^ t) / l;
        return PowerMod(g, pow, N);
    }

    ZZ hashprime(const ZZ &g, const ZZ &y) {
        std::string input = to_bin(g) + '*' + to_bin(y);
        return NextPrime(hash(input), NUM_TRIALS);
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

    RSA_Group group;

    ZZ x = RandomBits_ZZ(NUM_LENGTH);
    ZZ g = group.hash(to_bin(x));
    uint64_t time = 200;

    auto t1 = std::chrono::high_resolution_clock::now();

    std::cout << group.trapdoor(g, time) << std::endl << std::endl;

    auto t2 = std::chrono::high_resolution_clock::now();
    auto ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << ms.count() << std::endl;

    t1 = std::chrono::high_resolution_clock::now();

    ZZ s = group.eval(g, time);

    std::cout << s << std::endl << std::endl;

    t2 = std::chrono::high_resolution_clock::now();
    ms = std::chrono::duration_cast<std::chrono::milliseconds>(t2 - t1);
    std::cout << ms.count() << std::endl;

    t1 = std::chrono::high_resolution_clock::now();

    ZZ l = RandomPrime_ZZ(2 * NUM_LENGTH, 50);
    l = group.hashprime(g, s);

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
