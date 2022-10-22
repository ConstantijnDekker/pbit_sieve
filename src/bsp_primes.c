#include <stdint.h>
#include <bsp.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

uint64_t WHEEL_PRIMES[4] = { 2, 3, 5, 7 };
uint64_t WHEEL[48] = {
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181,
    187, 191, 193, 197, 199, 209,
};
uint64_t WHEEL_SZ = 210;
uint64_t BLOCK_SZ = 1 << 12;

// Plain sieve of Eratosthenes that finds all primes <= n that are not in the WHEEL.
uint64_t *simple_sieve(uint64_t n, uint64_t *nprimes) {
    bool *is_prime = malloc((n + 1) * sizeof(bool));
    for (uint64_t i = 0; i <= n; ++i) {
        is_prime[i] = true;
    }

    *nprimes = 0;
    for (uint64_t i = 2; i <= n; ++i) {
        if (is_prime[i]) {
            ++(*nprimes);
            for (uint64_t j = i * i; j <= n; j += i) {
                is_prime[j] = false;
            }
        }
    }
    uint64_t *sieving_primes = malloc(*nprimes * sizeof(uint64_t));
    uint64_t *ptr = sieving_primes;

    for (uint64_t i = 2; i <= n; ++i) {
        if (i == 2 || i == 3 || i == 5 || i == 7) {
            *nprimes -= 1;
            continue;
        }
        if (is_prime[i]) {
            *ptr = i;
            ++ptr;
        }
    }

    return sieving_primes;
}

// Set the flags that are k mod p in block. Return next offset.
uint64_t sieve_with_prime(uint64_t p, uint64_t k, uint64_t *block) {
    uint64_t j = k;
    while ((j >> 6) < BLOCK_SZ) {
        block[j >> 6] |= 1UL << (j & 0x3F);
        j += p;
    }
    return j - 64 * BLOCK_SZ;
}

// Sieve primes from a block.
void sieve_block(uint64_t *primes, uint64_t *offsets, uint64_t nprimes, uint64_t *block) {
    memset(block, 0, BLOCK_SZ * sizeof(uint64_t));
    for (uint64_t i = 0; i < nprimes; ++i) {
        offsets[i] = sieve_with_prime(primes[i], offsets[i], block);
    }
}

// Count the zero bits in a block.
uint64_t count_zeros(uint64_t *block) {
    uint64_t nzeros = BLOCK_SZ * 64;
    for (uint64_t i = 0; i < BLOCK_SZ; ++i) {
        nzeros -= (uint64_t) __builtin_popcountl(block[i]);
    }
    return nzeros;
}

uint64_t sieve_blocks(uint64_t nblocks, uint64_t *primes, uint64_t *offsets, uint64_t nprimes, uint64_t *block) {
    uint64_t prime_cnt = 0;
    while (nblocks--) {
        sieve_block(primes, offsets, nprimes, block);
        prime_cnt += count_zeros(block);
    }
    return prime_cnt;
}

// Compute a / b mod WHEEL_SZ
uint64_t mod_div(uint64_t a, uint64_t b) {
    for (uint64_t i = 0; i < WHEEL_SZ; ++i) {
        if ((WHEEL[i] * b) % WHEEL_SZ == a) {
            return WHEEL[i];
        }
    }
    // Error
    return 0;
}

// Count the number of primes congruent to w modulo WHEEL_SZ.
uint64_t count_res_class(uint64_t *primes, uint64_t nprimes, uint64_t nblocks, uint64_t w) {
    uint64_t block[BLOCK_SZ];
    uint64_t *offsets = malloc(nprimes * sizeof(uint64_t));
    for (uint64_t i = 0; i < nprimes; ++i) {
        offsets[i] = (mod_div(w, primes[i]) * primes[i]) / WHEEL_SZ;
    }
    return sieve_blocks(nblocks, primes, offsets, nprimes, block);
}

uint64_t count_primes(uint64_t nblocks, uint64_t *nprimes) {
    uint64_t n = nblocks * BLOCK_SZ * 64 * WHEEL_SZ;
    *nprimes = 0;
    uint64_t *primes = simple_sieve((uint64_t) sqrt((double) n), nprimes);
    uint64_t s = bsp_pid();

    uint64_t prime_cnt = 0;
    for (uint64_t i = 4 * s; i < 4 * (s + 1); ++i) {
        prime_cnt += count_res_class(primes, *nprimes, nblocks, WHEEL[i]);
    }

    // Add the primes in the wheel that were not tested.
    // Also add primes used for sieving because they sieved themselves out.
    // Finally subtract 1 because it isn't prime.
    return prime_cnt;
}

unsigned P = 12;

void bsp_prime_count() {
    bsp_begin(P);
    uint64_t nblocks = 500;
    uint64_t nprimes = 0;
    uint64_t prime_cnt = count_primes(nblocks, &nprimes);
    bsp_pid_t s = bsp_pid();
    uint64_t counts[P];
    bsp_push_reg(counts, P * sizeof(uint64_t));
    bsp_sync();
    for (bsp_pid_t t = 0; t < P; ++t) {
        bsp_put(t, &prime_cnt, counts, s * sizeof(uint64_t), sizeof(uint64_t));
    }
    bsp_sync();
    uint64_t total = 0;
    for (bsp_pid_t t = 0; t < P; ++t) {
        total += counts[t];
    }
    total = total + 4 + nprimes - 1;

    if (s == 0) {
        printf("%lu primes below %lu\n", total, nblocks * BLOCK_SZ * 64 * WHEEL_SZ);
    }
}

int main(int argc, char *argv[]) {
    bsp_init(&bsp_prime_count, argc, argv);
    bsp_prime_count();
    return 0;
}

