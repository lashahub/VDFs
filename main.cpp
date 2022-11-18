#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include "NTL/ZZ.h"
#include "bicycl.hpp"
#include "openssl/evp.h"
#include "openssl/sha.h"

#define NUM_TRIALS 50
#define SECURITY_PARAM 128
#define NUM_LENGTH 1536

using namespace NTL;
using namespace BICYCL;

class RSA_Group {

public:
    ZZ p, q, N, fi;

public:

    RSA_Group() {
        RandomPrime(this->p, NUM_LENGTH, NUM_TRIALS);
        RandomPrime(this->q, NUM_LENGTH, NUM_TRIALS);
        N = p * q;
        fi = (p - 1) * (q - 1);
    }

    ZZ bytes_to_ZZ(std::vector<uint8_t> bytes) {
        std::reverse(bytes.begin(), bytes.end());
        return ZZFromBytes(&bytes[0], bytes.size());
    }

    std::vector<uint8_t> ZZ_to_bytes(const ZZ &x) {
        long len = NumBytes(x);
        std::vector<uint8_t> bytes(len);
        BytesFromZZ(&bytes[0], x, len);
        std::reverse(bytes.begin(), bytes.end());
        return bytes;
    }

    std::string bytes_to_hex_string(const std::vector<uint8_t> &bytes) {
        std::ostringstream stream;
        for (uint8_t b: bytes) {
            stream << std::setw(2) << std::setfill('0') << std::hex << static_cast<int>(b);
        }
        return stream.str();
    }

    std::vector<uint8_t> SHA256(const std::string &input) {
        uint32_t digest_length = SHA256_DIGEST_LENGTH;
        const EVP_MD *algorithm = EVP_sha3_256();
        auto *digest = static_cast<uint8_t *>(OPENSSL_malloc(digest_length));
        EVP_MD_CTX *context = EVP_MD_CTX_new();
        EVP_DigestInit_ex(context, algorithm, nullptr);
        EVP_DigestUpdate(context, input.c_str(), input.size());
        EVP_DigestFinal_ex(context, digest, &digest_length);
        EVP_MD_CTX_destroy(context);
        std::vector<uint8_t> bytes(digest, digest + digest_length);
        OPENSSL_free(digest);
        return bytes;
    }

    ZZ hash(const ZZ &x) {
        std::string input = "residue" + bytes_to_hex_string(ZZ_to_bytes(x));
        return bytes_to_ZZ(SHA256(input)) % N;
    }


    ZZ trapdoor(const ZZ &g, long t) const {
        ZZ e = PowerMod(ZZ(2), t, fi);
        return PowerMod(g, e, N);
    }

    ZZ eval(const ZZ &g, long t) const {
        ZZ e = power2_ZZ(t);
        return PowerMod(g, e, N);
    }

    ZZ random_prime() {
        return RandomPrime_ZZ(2 * (SECURITY_PARAM + 10), NUM_TRIALS);
    }

    ZZ noninteractive_prime(const ZZ &g, const ZZ &y) {
        auto g_s = ZZ_to_bytes(g);
        g_s.push_back('*');
        auto y_s = ZZ_to_bytes(y);
        g_s.insert(g_s.end(), y_s.begin(), y_s.end());
        ZZ numeric = bytes_to_ZZ(g_s);
        numeric = numeric % power2_ZZ(2 * (SECURITY_PARAM + 10));
        return NextPrime(numeric, NUM_TRIALS);
    }

    ZZ trapdoor_proof(const ZZ &g, long t, const ZZ &l) {
        ZZ e = (power2_ZZ(t) / l) % fi;
        return PowerMod(g, e, N);
    }

    ZZ eval_proof(const ZZ &g, long t, const ZZ &l) const {
        ZZ x(1);
        ZZ r(1);
        for (long i = 0; i < t; i++) {

        }
        return x;
    }

    bool verify(const ZZ &g, const ZZ &pi, const ZZ &l, const ZZ &r, const ZZ &y) {
        ZZ fst = PowerMod(pi, l, N);
        ZZ snd = PowerMod(g, r, N);
        return MulMod(fst, snd, N) == y;
    }

};

class Class_Group {

public:
    ClassGroup *cg;
    Mpz disc;

public:

    Class_Group() {
        RandGen r;
        disc = r.random_negative_discriminant(NUM_LENGTH);
        cg = new ClassGroup(disc);
    }

    ~Class_Group() {
        delete cg;
    }

    QFI eval(const QFI &g, long t) {
        Mpz e;
        mpz_ui_pow_ui(e.mpz_, 2, t);
        QFI res;
        cg->nupow(res, g, e);
        return res;
    }


};

int main() {
    Class_Group cg;
    RandGen r;
    QFI g = cg.cg->random(r);

    while (true) {

        long n, b;
        std::cin >> n;

        auto t01 = std::chrono::high_resolution_clock::now();
        std::cout << cg.eval(g, n) << std::endl;
        auto t02 = std::chrono::high_resolution_clock::now();
        auto m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
        std::cout << m0s.count() << std::endl;

    }

}

int main1() {

    RSA_Group group;

    ZZ x = RandomBits_ZZ(NUM_LENGTH);
    ZZ g = group.hash(x);
    long time = 2000000;

    ZZ y = group.trapdoor(g, time);
    std::cout << y << std::endl << std::endl;

    auto t01 = std::chrono::high_resolution_clock::now();
    std::cout << group.eval(g, time) << std::endl << std::endl;
    auto t02 = std::chrono::high_resolution_clock::now();
    auto m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << m0s.count() << std::endl;

    ZZ l = group.random_prime();

    ZZ pi = group.trapdoor_proof(g, time, l);
    ZZ r = PowerMod(ZZ(2), time, l);
    std::cout << group.verify(g, pi, l, r, y) << std::endl;

    return 0;
}
