#include <iostream>
#include <cmath>
#include <chrono>
#include <iomanip>
#include <omp.h>
#include "NTL/ZZ.h"
#include "bicycl.hpp"
#include "openssl/evp.h"
#include "openssl/sha.h"

#define NUM_TRIALS 50
#define SECURITY_PARAM 128
#define NUM_LENGTH 1536
#define NUM_THREADS 12

using namespace NTL;
using namespace BICYCL;

class RSA_Group {

public:
    ZZ p, q, N, fi, lambda;

public:

    RSA_Group() {
        RandomPrime(this->p, NUM_LENGTH, NUM_TRIALS);
        RandomPrime(this->q, NUM_LENGTH, NUM_TRIALS);
        N = p * q;
        fi = (p - 1) * (q - 1);
        fi /= GCD(p - 1, q - 1);
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
        return bytes_to_ZZ(SHA256(input));
    }

    ZZ noninteractive_prime(const ZZ &g, const ZZ &y) {
        auto g_s = ZZ_to_bytes(g);
        g_s.push_back('*');
        auto y_s = ZZ_to_bytes(y);
        g_s.insert(g_s.end(), y_s.begin(), y_s.end());
        ZZ numeric = bytes_to_ZZ(g_s);
        numeric = hash(numeric);
        //numeric = numeric % power2_ZZ(2 * (SECURITY_PARAM + 10));
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
        ZZ pwr = ((1 << K) * (power2_ZZ(t - K * (i + 1)) % l)) / l;
        long res;
        conv(res, pwr);
        return res;
    }

    ZZ fast_mul(std::vector<ZZ> &list, std::vector<long> &indices) {
        ZZ cores[NUM_THREADS];
        for (auto &core: cores) {
            core = ZZ(1);
        }

#pragma omp parallel for
        for (long i = 0; i < indices.size(); i++) {
            int n = omp_get_thread_num();
            cores[n] = MulMod(cores[n], list[indices[i]], N);
        }
        ZZ res(1);
        for (long i = 0; i < NUM_THREADS; i++) {
            res = MulMod(res, cores[i], N);
        }
        return res;
    }

    // Algorithm 3 & 5.2
    std::pair<ZZ, ZZ> eval_fast(const ZZ &x, long t) {

        long K = (long) log2((double) t) / 2;
        ZZ g = hash(x);
        std::vector<ZZ> C;
        ZZ y = g;

        for (long i = 0; i < t; i++) {
            if (i % K == 0) C.push_back(y);
            y = SqrMod(y, N);
        }

        ZZ l = noninteractive_prime(g, y);

        omp_set_num_threads(NUM_THREADS);
        std::vector<std::vector<long>> Is_T[NUM_THREADS];
#pragma omp parallel for
        for (long i = 0; i < NUM_THREADS; i++) {
            Is_T[i] = std::vector<std::vector<long>>(1 << K);
            for (long j = 0; j < (1 << K); j++) {
                Is_T[i][j] = std::vector<long>();
            }
        }

#pragma omp parallel for
        for (long i = 0; i < t / K; i++) {
            long b_i = get_block(t, K, i, l);
            Is_T[omp_get_thread_num()][b_i].push_back(i);
        }

        std::vector<std::vector<long>> Is(1 << K);

#pragma omp parallel for
        for (long i = 0; i < (1 << K); i++) {
            for (auto &j: Is_T) {
                Is[i].reserve(Is[i].size() + distance(j[i].begin(), j[i].end()));
                Is[i].insert(Is[i].end(), j[i].begin(), j[i].end());
            }
        }

        ZZ pi(1);
        for (long b = 0; b < (1 << K); b++) {
            auto Ib = Is[b];
            ZZ prod = fast_mul(C, Ib);
            ZZ pwr = PowerMod(prod, b, N);
            pi = MulMod(pi, pwr, N);
        }

        return std::make_pair(y, pi);

    }

    // Algorithm 3 & 5.1
    std::pair<ZZ, ZZ> eval_fast_highmem(const ZZ &x, long t) {
        long K = (long) log2((double) t) / 2;
        std::cout << K << std::endl;
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

    ZZ random_prime() {
        return RandomPrime_ZZ(264, NUM_TRIALS);
    }

    std::pair<QFI, QFI> eval(const QFI &g, long t, const ZZ &l) {

        auto t01 = std::chrono::high_resolution_clock::now();
        Mpz e;
        mpz_ui_pow_ui(e.mpz_, 2, t);
        QFI y;
        cg->nupow(y, g, e);
        auto t02 = std::chrono::high_resolution_clock::now();
        auto m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
        std::cout << m0s.count() << std::endl;

        ZZ e1 = power2_ZZ(t) / l;
        Mpz em = ZZ_to_mpz(e1);
        QFI pi;
        cg->nupow(pi, g, em);
        return std::make_pair(y, pi);
    }

    Mpz ZZ_to_mpz(const ZZ &x) {
        Mpz res;
        std::stringstream ss1;
        ss1 << x;
        mpz_set_str(res.mpz_, ss1.str().c_str(), 10);
        return res;
    }

    bool verify(const QFI &x, const QFI &y, const QFI &pi, long t, const ZZ &l) {
        ZZ r = power2_ZZ(t) % l;

        QFI fst;
        cg->nupow(fst, pi, ZZ_to_mpz(l));
        QFI snd;
        cg->nupow(snd, x, ZZ_to_mpz(r));
        QFI res;
        cg->nucomp(res, fst, snd);
        return res == y;
    }


};

int main_class_groups() {
    Class_Group cg;
    RandGen r;
    QFI g = cg.cg->random(r);

    long time = 2000000;

    auto t01 = std::chrono::high_resolution_clock::now();
    ZZ l = cg.random_prime();
    auto res = cg.eval(g, time, l);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << m0s.count() << std::endl;
    std::cout << ((cg.verify(g, res.first, res.second, time, l)) ? "Verified" : "Failed") << std::endl;

    return 0;

}

int main() {

    RSA_Group group;

    ZZ x = RandomBits_ZZ(2 * NUM_LENGTH);
    long time = 2000000;

    //////////////////////////////////////////////////////////////////////////////

    auto t01 = std::chrono::high_resolution_clock::now();
    auto res = group.trapdoor(x, time);
    auto t02 = std::chrono::high_resolution_clock::now();
    auto m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << m0s.count() << std::endl;
    std::cout << res.first << std::endl;
    std::cout << res.second << std::endl;
    std::cout << (group.verify(x, res.first, res.second, time) ? "Verified" : "Failed") << std::endl;

    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////

    t01 = std::chrono::high_resolution_clock::now();
    res = group.eval(x, time);
    t02 = std::chrono::high_resolution_clock::now();
    m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << m0s.count() << std::endl << std::endl;
    std::cout << res.first << std::endl;
    std::cout << res.second << std::endl;
    std::cout << (group.verify(x, res.first, res.second, time) ? "Verified" : "Failed") << std::endl;

    //////////////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////////////

    t01 = std::chrono::high_resolution_clock::now();
    res = group.eval_fast(x, time);
    t02 = std::chrono::high_resolution_clock::now();
    m0s = std::chrono::duration_cast<std::chrono::milliseconds>(t02 - t01);
    std::cout << m0s.count() << std::endl;
    std::cout << res.first << std::endl;
    std::cout << res.second << std::endl;
    std::cout << (group.verify(x, res.first, res.second, time) ? "Verified" : "Failed") << std::endl;

    //////////////////////////////////////////////////////////////////////////////

    return 0;

}
