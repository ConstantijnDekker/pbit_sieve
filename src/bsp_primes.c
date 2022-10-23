#include <stdint.h>
#include <bsp.h>
#include <stdio.h>
#include <stdbool.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

const uint64_t WHEEL_PRIMES[4] = { 2, 3, 5, 7};
const size_t WHEEL_LEN = 48;
const uint64_t WHEEL[48] = {
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181,
    187, 191, 193, 197, 199, 209,
};
const size_t WORD_SZ = 64;
const size_t WHEEL_SZ = 210;
const size_t BLOCK_SZ = 1 << 12;

// Plain sieve of Eratosthenes that finds all primes <= n that are not in the WHEEL.
uint64_t *simple_sieve(size_t n, size_t *nprimes) {
    bool *is_prime = malloc((n + 1) * sizeof(bool));
    for (size_t i = 0; i <= n; ++i) {
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
    free(is_prime);

    return sieving_primes;
}

// Set the flags that are k mod p in block. Return next offset.
size_t sieve_with_prime(uint64_t p, size_t k, size_t *block, size_t block_len) {
    size_t j = k;
    while (j < WORD_SZ * block_len) {
        block[j >> 6] |= 1UL << (j & 0x3F);
        j += p;
    }
    return j - WORD_SZ * block_len;
}

// Sieve primes from a block.
void sieve_block(uint64_t *primes, size_t *offsets, size_t nprimes, size_t *block, size_t block_len) {
    memset(block, 0, BLOCK_SZ * sizeof(size_t));
    for (size_t i = 0; i < nprimes; ++i) {
        offsets[i] = sieve_with_prime(primes[i], offsets[i], block, block_len);
    }
}

uint64_t zero_cnt(uint64_t word) {
    return WORD_SZ - __builtin_popcountl(word);
}

// Count the zero bits in a block.
uint64_t count_zeros(size_t *block, size_t block_len) {
    uint64_t nzeros = 0;
    for (size_t i = 0; i < block_len; ++i) {
        nzeros += zero_cnt(block[i]);
    }
    return nzeros;
}

uint64_t mask(uint64_t width) {
    return (1ULL << width) - 1;
}

uint64_t sieve_blocks(size_t nbits, uint64_t *primes, size_t *offsets, size_t nprimes, size_t *block) {
    uint64_t prime_count = 0;
    while (nbits >= WORD_SZ * BLOCK_SZ) {
        sieve_block(primes, offsets, nprimes, block, BLOCK_SZ);
        prime_count += count_zeros(block, BLOCK_SZ);
        nbits -= WORD_SZ * BLOCK_SZ;
    }
    
    // The final block does not have to be fully processed.
    size_t last = nbits / WORD_SZ;
    sieve_block(primes, offsets, nprimes, block, last + 1);
    prime_count += count_zeros(block, last + 1);

    // small correction for overcounting zeros in the last word
    return prime_count - zero_cnt(block[last] | mask(nbits % WORD_SZ));
}

// Compute a / b mod WHEEL_SZ
uint64_t mod_div(uint64_t a, uint64_t b) {
    for (size_t i = 0; i < WHEEL_SZ; ++i) {
        if ((WHEEL[i] * b) % WHEEL_SZ == a) {
            return WHEEL[i];
        }
    }
    // Error
    return 0;
}

// Count the number of primes congruent to w modulo WHEEL_SZ.
uint64_t count_res_class(size_t nbits, uint64_t *primes, size_t nprimes, uint64_t w) {
    size_t block[BLOCK_SZ];
    size_t *offsets = malloc(nprimes * sizeof(size_t));
    for (size_t i = 0; i < nprimes; ++i) {
        uint64_t p = primes[i];
        //offsets[i] = (mod_div(w, primes[i]) * primes[i]) / WHEEL_SZ;
        uint64_t r = mod_div(w, p);
        // Make sure r >= p while remaining the same mod WHEEL_SZ
        if (r >= p % WHEEL_SZ) {
            r += p - (p % WHEEL_SZ);
        } else {
            r += p - (p % WHEEL_SZ) + WHEEL_SZ;
        }
        // bit index of first multiple of p to be crossed of
        offsets[i] = (r * p) / WHEEL_SZ;
    }
    uint64_t prime_count = sieve_blocks(nbits, primes, offsets, nprimes, block);
    free(offsets);
    return prime_count;
}

uint64_t count_primes(size_t n) {
    uint64_t nprimes = 0;
    uint64_t *primes = simple_sieve((uint64_t) sqrt((double) n), &nprimes);
    bsp_pid_t s = bsp_pid();
    bsp_pid_t p = bsp_nprocs();

    uint64_t prime_count = 0;
    for (size_t i = s; i < WHEEL_LEN; i += p) {
        size_t nbits = n / WHEEL_SZ;
        if (WHEEL[i] <= n % WHEEL_SZ) {
            nbits += 1;
        }
        prime_count += count_res_class(nbits, primes, nprimes, WHEEL[i]);
    }

    free(primes);
    return prime_count;
}

uint64_t P = 12;
uint64_t N = 0;

void bsp_prime_count() {
    bsp_begin(P);
    size_t n = N;
    uint64_t prime_count = count_primes(n);
    bsp_pid_t s = bsp_pid();
    uint64_t counts[P];
    bsp_push_reg(counts, P * sizeof(uint64_t));
    bsp_sync();
    for (bsp_pid_t t = 0; t < P; ++t) {
        bsp_put(t, &prime_count, counts, s * sizeof(uint64_t), sizeof(uint64_t));
    }
    bsp_sync();
    uint64_t total = 0;
    for (bsp_pid_t t = 0; t < P; ++t) {
        total += counts[t];
    }

    for (size_t i = 0; i < 4; ++i) {
        if (WHEEL_PRIMES[i] <= n) {
            total += 1;
        }
    }

    if (n >= 1) {
        total -= 1;
    }

    if (s == 0) {
        printf("%lu primes below %lu\n", total, n);
    }
    bsp_end();
}

int main(int argc, char *argv[]) {
    if (argc < 2 || argc > 3) {
        printf("Incorrect amount of arguments\n");
        return 0;
    }
    bsp_init(&bsp_prime_count, argc, argv);

    sscanf(argv[1], "%lu", &N);
    if (argc == 3) {
        sscanf(argv[2], "%lu", &P);
    }
    bsp_prime_count();
    return 0;
}

