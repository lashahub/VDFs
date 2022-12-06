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

    // Algorithm 1
    std::pair<ZZ, ZZ> trapdoor(const ZZ &x, long t) {
        ZZ g = hash(x);
        ZZ e = PowerMod(ZZ(2), t, fi);
        ZZ y = PowerMod(g, e, N);

        ZZ l = noninteractive_prime(g, y);
        ZZ rq = (power2_ZZ(t) / l) % fi;
        ZZ pi = PowerMod(g, rq, N);

        return std::make_pair(y, pi);
    }

    //Algorithm 2
    bool verify(const ZZ &x, const ZZ &y, const ZZ &pi, long t) {
        ZZ g = hash(x);
        ZZ l = noninteractive_prime(g, y);
        ZZ r = power2_ZZ(t) % l;

        ZZ fst = PowerMod(pi, l, N);
        ZZ snd = PowerMod(g, r, N);
        return MulMod(fst, snd, N) == y;
    }

    // Algorithm 3
    std::pair<ZZ, ZZ> eval(const ZZ &x, long t) {
        ZZ g = hash(x);
        ZZ y = PowerMod(g, power2_ZZ(t), N);

        ZZ l = noninteractive_prime(g, y);
        ZZ pi = simple_proof(g, l, t);
        return std::make_pair(y, pi);
    }

    // Algorithm 4
    ZZ simple_proof(const ZZ &g, const ZZ &l, long t) {
        ZZ e = power2_ZZ(t) / l;
        return PowerMod(g, e, N);
    }

    // Algorithm 4 Paper
    ZZ simple_proof_paper(const ZZ &g, const ZZ &l, long t) {
        ZZ x(1);
        ZZ r(1);
        for (long i = 0; i < t; i++) {
            ZZ b = 2 * r / l;
            r = 2 * r % l;
            x = MulMod(SqrMod(x, N), PowerMod(g, b, N), N);
        }
        return x;
    }

    long get_block(long t, long K, long i, const ZZ &l) {
        ZZ pwr = (power2_ZZ(K) * (power2_ZZ(t - K * (i + 1)) % l)) / l;
        long res;
        conv(res, pwr);
        return res;
    }

    // Algorithm 3 & 5.1
    std::pair<ZZ, ZZ> eval_fast_highmem(const ZZ &x, long t) {
        long K = (long) log2((double) t) / 2;
        std::cout<<K<<std::endl;
        ZZ g = hash(x);
        std::vector<ZZ> C;

        ZZ y = g;

        auto t01 = std::chrono::high_resolution_clock::now();

        for (long i = 0; i < t; i++) {
            if (i % K == 0) C.push_back(y);
            y = SqrMod(y, N);
        }

        auto t02 = std::chrono::high_resolution_clock::now();
        auto m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
        std::cout << m0s.count() << std::endl;

        ZZ l = noninteractive_prime(g, y);

        ZZ total_power = power2_ZZ(t) / l;


        t01 = std::chrono::high_resolution_clock::now();

        ZZ pi(1);
        for (long i = C.size() - 1; i >= 0; i--) {
            //if (i % 10000 == 0) std::cout << i << std::endl;
            ZZ curr_power = power2_ZZ(K * i);
            pi = MulMod(pi, PowerMod(C[i], total_power / curr_power, N), N);
            total_power %= curr_power;
        }


        t02 = std::chrono::high_resolution_clock::now();
        m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
        std::cout << m0s.count() << std::endl;

        return std::make_pair(y, pi);

    }

    // Algorithm 3 & 5
    std::pair<ZZ, ZZ> eval_fast_lowmem(const ZZ &x_orig, long t) {

        long K = (long) log2((double) t) / 3;
        long gamma = (long) sqrt((double) t);
        long Kg = K * gamma;
        long ceiling = t / Kg + (t % Kg != 0);

        ZZ g = hash(x_orig);
        std::vector<ZZ> C(ceiling + 1);
        ZZ tmp = g;

        ZZ y;

        for (long i = 0; i <= t + Kg; i++) {
            if (i == t) y = tmp;
            if (i % Kg == 0) {
                C[i / Kg] = tmp;
            }
            tmp = SqrMod(tmp, N);
        }

        ZZ l = noninteractive_prime(g, y);

        long K1 = K / 2;
        long K0 = K - K1;
        ZZ x(1);
        for (long j = gamma - 1; j >= 0; j--) {
            x = PowerMod(x, power2_ZZ(K), N);
            std::vector<ZZ> Y(1 << K);
            for (long b = 0; b < (1 << K); b++) {
                Y[b] = ZZ(1);
            }

            for (long i = 0; i < ceiling; i++) {
                long b = get_block(t, K, i * gamma + j, l);
                Y[b] = MulMod(Y[b], C[i], N);
            }

            for (long b1 = 0; b1 < (1 << K1); b1++) {
                ZZ z = ZZ(1);
                for (long b0 = 0; b0 < (1 << K0); b0++) {
                    z = MulMod(z, Y[b1 * (1 << K0) + b0], N);
                }
                x = MulMod(x, PowerMod(z, b1 * (1 << K0), N), N);
            }

            for (long b0 = 0; b0 < (1 << K0); b0++) {
                ZZ z = ZZ(1);
                for (long b1 = 0; b1 < (1 << K1); b1++) {
                    z = MulMod(z, Y[b1 * (1 << K0) + b0], N);
                }
                x = MulMod(x, PowerMod(z, b0, N), N);
            }

        }
        return std::make_pair(y, x);
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

int main1() {
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

int main() {

    RSA_Group group;

    ZZ x = RandomBits_ZZ(NUM_LENGTH);
    long time = 2000000;

    //////////////////////////////////////////////////////////////////////////////

    auto t01 = std::chrono::high_resolution_clock::now();
    group.trapdoor(x, time);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << m0s.count() << std::endl;

    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////

    t01 = std::chrono::high_resolution_clock::now();
    group.eval(x, time);
    t02 = std::chrono::high_resolution_clock::now();
    m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << m0s.count() << std::endl << std::endl;

    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////

    t01 = std::chrono::high_resolution_clock::now();
    group.eval_fast_highmem(x, time);
    t02 = std::chrono::high_resolution_clock::now();
    m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << m0s.count() << std::endl;

    //////////////////////////////////////////////////////////////////////////////

}
