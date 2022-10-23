// Version from March 2021
// compiling with: g++ -Ofast -funroll-all-loops -ftree-vectorize -march=native ./filename
// yields 160 ms for primes below 10^9
#include <vector>
#include <iostream>
#include <algorithm>
#include <cstring>
#include <cmath>
#include <cassert>


// Print a vector with spaces between.
// Useful when debugging.
#define printVector(v) for (int i = 0; i < (int)(v).size(); ++i) { cout << (v)[i] << ' '; } cout << '\n';

#define BIT0 1
#define BIT1 2
#define BIT2 4
#define BIT3 8
#define BIT4 16
#define BIT5 32
#define BIT6 64
#define BIT7 128

// Saves typing.
using namespace std;

/* Program constants. */
const int block_size = (1 << 15);
const vector<int64_t> wheel = {1, 7, 11, 13, 17, 19, 23, 29};
const vector<int64_t> gaps  = {2, 6,  4,  2,  4,  2,  4,  6};
/* wheel[wheel_loc[i]] == highest number in the wheel
 * that is less than or equal to i, for 0 <= i < 30. */
const vector<int64_t> wheel_loc = {0, 0, 0, 0, 0, 0, 0, 1, 1, 1, 1, 2, 2,
    3, 3, 3, 3, 4, 4, 5, 5, 5, 5, 6, 6, 6, 6, 6, 6, 7};
const int64_t wheelsz = 8;
const int64_t Pi = 30;
// product of presieve primes.
const int64_t preS = 7 * 11 * 13 * 17 * 19;
// number of primes we pre_sieve.
const int64_t pre_sieve_nr = 5;

// Copy content from the presieve into our block.
void init_block(uint8_t* block, uint8_t* pre_sieve, int64_t& loc) {
    // index in block to copy to.
    int64_t q = 0;
    while (q < block_size) {
        int64_t x = min(preS - loc, block_size - q);
        // preS == loc, so loc as to start back at beginning.
        if (x == 0) {
            loc = 0;
            continue;
        }
        memcpy(block + q, pre_sieve + loc, x * sizeof(uint8_t));
        loc += x;
        q += x;
    }
}

/* Simple sieve to find all primes below or equal to S.
 * The primes will be used for a better sieve. */
vector<int64_t> get_sieving_primes(int64_t S) {
    vector<int64_t> sieving_primes;
    int64_t lS = (int64_t) log(S);
    sieving_primes.reserve(S / lS);
    vector<char> is_composite(S + 1, 0);
    for (int64_t i = 3; i <= S; i += 2) {
        if (not is_composite[i]) {
            /* Do not add the primes 3 and 5,
             * we will not need them. */
            if (i != 3 && i != 5) {
                sieving_primes.push_back(i);
            }
            for (int64_t j = i * i; j <= S; j += i) {
                is_composite[j] = 1;
            }
        }
    }
    return sieving_primes;
}

// For every sieving prime. Store the byte index
// of the start of the wheel where it will start marking things off.
void set_next_bytes(const vector<int64_t>& primes, const vector<int64_t>& wheel_start,
        const vector< vector<int64_t> >& byte_sums, vector<int64_t>& next_byte) {
    for (int64_t i = 0; i < (int64_t) primes.size(); ++i) {
        int64_t S2 = primes[i] * primes[i];
        next_byte[i] = S2 / Pi;
        next_byte[i] -= byte_sums[i][wheel_start[i]];
    }
}

/* For every sieving prime. Store where it starts
 * in our wheel. This is to avoid divisions later. */
void set_wheel_index(const vector<int64_t>& primes, vector<int64_t>& wheel_index) {
    for (int64_t i = 0; i < (int64_t) primes.size(); ++i) {
        wheel_index[i] = wheel_loc[primes[i] % Pi];
    }
}

// pattern[i][j] is defined to
// be a bit in the position of wheel[i] * wheel[j].
void set_pattern(uint8_t pattern[8][8]) {
    for (int64_t i = 0; i < 8; ++i) {
        for (int64_t j = 0; j < wheelsz; ++j) {
            pattern[i][j] = 1 << wheel_loc[((wheel[i] * wheel[j]) % Pi)];
        }
    }
}

// Calculate gaps between markings in bytes.
// Really complicated stuff.
void set_byte_gaps(const vector<int64_t>& primes, uint8_t pattern[8][8],
        const vector<int64_t>& wheel_start, vector< vector<int64_t> >& byte_gaps) {
    for (int64_t i = 0; i < (int64_t) primes.size(); ++i) {
        int64_t div = primes[i] / Pi;
        for (int64_t j = 1; j < 8; ++j) {
            byte_gaps[i][j - 1] = div * gaps[j] + ((primes[i] % Pi) * gaps[j]) / Pi;
            byte_gaps[i][j - 1] += 
                pattern[wheel_start[i]][j] < pattern[wheel_start[i]][j - 1] ? 1 : 0;
        }
        byte_gaps[i][7] = div * gaps[0] + ((primes[i] % Pi) * gaps[0]) / Pi;
        byte_gaps[i][7] += 
            pattern[wheel_start[i]][0] < pattern[wheel_start[i]][wheelsz - 1] ? 1 : 0;
    }
}

// Calculate an offset for each new mark to be made.
void set_byte_sums(const vector< vector<int64_t> >& byte_gaps,
        vector<vector<int64_t>>& byte_sums) {
    for (int64_t i = 0; i < (int64_t) byte_gaps.size(); ++i) {
        byte_sums[i][0] = 0;
        for (int64_t j = 1; j < wheelsz; ++j) {
            byte_sums[i][j] = byte_sums[i][j - 1] + byte_gaps[i][j - 1];
        }
    }
}

// Convenience lookup table. for 1 <= x < 256, unset_bits[x] is the
// number of zeroes in the binary expansion of x.
void build_zeroes_table(vector<int64_t>& unset_bits) {
    unset_bits[0] = 8;
    for (int64_t i = 0; i < unset_bits.size(); ++i) {
        unset_bits[i] = unset_bits[i / 2] - (i & 1);
    }
}

// Function that adds up the number of bits set to zero in the block.
// Warning: only use when block_size is divisible by 8!
int64_t count_primes(uint8_t block[]) {
    int64_t prime_count = 0;
    for (int64_t i = 0; i < block_size; i += 8) {
        prime_count += __builtin_popcountll(~(*(reinterpret_cast<uint64_t*>(block + i))));
    }
    return prime_count;
}


int main(int argc, char *argv[]) {
    long long N;
    // By default, run the program on an input size
    // of 1 billion :). It can handle this input in 160ms!.
    if (argc <= 1) {
        N = 1'000'000'000;
    } else {
        string s(argv[1]);
        N = stoll(s);
    }
    if (N < 0) {
        return 1;
    }

    int64_t S = (int64_t) sqrtl(N) + 1;
    // Make sure we find at least the primes
    // up to 19. This is necessary to avoid segmentation fault.
    S = max(S, (int64_t) 400);

    /* Vector of primes we will use to sieve.
     * it starts at 7. */
    const vector<int64_t> primes = get_sieving_primes(S);
    const int64_t ss = primes.size();

    /* Allocate a block in which we will mark off multiples. */
    uint8_t block[block_size];

    vector<int64_t> squares(ss);
    for (int i = 0; i < ss; ++i) {
        squares[i] = (primes[i] * primes[i]) / Pi;
    }

    /* For every sieving prime, store the
     * index of where it sits in the wheel. */
    vector<int64_t> wheel_index(ss);
    set_wheel_index(primes, wheel_index);
    const vector<int64_t> wheel_start = wheel_index;

    uint8_t pattern[wheelsz][wheelsz];
    set_pattern(pattern);

    // Calculation of the gaps. Never look at it! :)
    vector< vector<int64_t> > byte_gaps(ss, vector<int64_t>(wheelsz));
    set_byte_gaps(primes, pattern, wheel_start, byte_gaps);

    // Summation over the gaps.
    vector< vector<int64_t> > byte_sums(ss, vector<int64_t>(wheelsz));
    set_byte_sums(byte_gaps, byte_sums);

    // store index of next_byte to mark off.
    vector<int64_t> next_byte(ss);
    set_next_bytes(primes, wheel_start, byte_sums, next_byte);

    vector<int64_t> unset_bits(256);
    build_zeroes_table(unset_bits);

    int64_t prime_count = 0;
    int64_t candidates = N / Pi * wheelsz + wheel_loc[N % Pi]
        + ((N % Pi) ? 1 : 0);
    int64_t final_index = candidates >> 3;

    // Calculate presieved pattern
    // based on the presieve primes.
    uint8_t* pre_sieve = (uint8_t*) calloc(preS, sizeof(uint8_t));
    for (int64_t i = 0; i < pre_sieve_nr; ++i) {
        int64_t indexx = primes[i] / Pi;
        int64_t k = 0;
        int64_t t = wheel_start[i];
        while (indexx < preS) {
            pre_sieve[indexx] |= pattern[t][k];
            indexx += byte_gaps[i][k];
            ++k; k &= 7;
        }
    }

    // Global index to indicate offset from where we began.
    int64_t offset = 0;
    // Index to keep track of where we are in the pre_sieve.
    int64_t ptr = 0;
    int64_t total = 0;
    while (offset <= final_index) {
        // Put the presieved pattern into the block.
        init_block(block, pre_sieve, ptr);
        // Set limit for which primes we want to search.
        int64_t lim = lower_bound(squares.begin(), squares.end(),
                offset + block_size) - squares.begin();
        int64_t index = 0;
        int64_t k = 0;
        for (int64_t i = pre_sieve_nr; i < lim; ++i) {
            int64_t bound = block_size - byte_sums[i][7];
            int64_t prime_type;
            index = next_byte[i] - offset;
            k = wheel_index[i];
            if (index < bound) {
                // Fall through until k is zero.
                switch(k) {
                    case 1:
                        block[index + byte_sums[i][1]] |= pattern[wheel_start[i]][1];
                    case 2:
                        block[index + byte_sums[i][2]] |= pattern[wheel_start[i]][2];
                    case 3:
                        block[index + byte_sums[i][3]] |= pattern[wheel_start[i]][3];
                    case 4:
                        block[index + byte_sums[i][4]] |= pattern[wheel_start[i]][4];
                    case 5:
                        block[index + byte_sums[i][5]] |= pattern[wheel_start[i]][5];
                    case 6:
                        block[index + byte_sums[i][6]] |= pattern[wheel_start[i]][6];
                    case 7:
                        block[index + byte_sums[i][7]] |= pattern[wheel_start[i]][7];
                        index += primes[i];
                }
                switch(wheel_start[i]) {
                    // p == 1 mod 30
                    case 0:
                        while (index < bound) {
                            block[index]                   |= BIT0;
                            block[index + byte_sums[i][1]] |= BIT1;
                            block[index + byte_sums[i][2]] |= BIT2;
                            block[index + byte_sums[i][3]] |= BIT3;
                            block[index + byte_sums[i][4]] |= BIT4;
                            block[index + byte_sums[i][5]] |= BIT5;
                            block[index + byte_sums[i][6]] |= BIT6;
                            block[index + byte_sums[i][7]] |= BIT7;
                            index += primes[i];
                        }
                        break;
                    // p == 7 mod 30
                    case 1:
                        while (index < bound) {
                            block[index]                   |= BIT1;
                            block[index + byte_sums[i][1]] |= BIT5;
                            block[index + byte_sums[i][2]] |= BIT4;
                            block[index + byte_sums[i][3]] |= BIT0;
                            block[index + byte_sums[i][4]] |= BIT7;
                            block[index + byte_sums[i][5]] |= BIT3;
                            block[index + byte_sums[i][6]] |= BIT2;
                            block[index + byte_sums[i][7]] |= BIT6;
                            index += primes[i];
                        }
                        break;
                    // p == 11 mod 30
                    case 2:
                        while (index < bound) {
                            block[index]                   |= BIT2;
                            block[index + byte_sums[i][1]] |= BIT4;
                            block[index + byte_sums[i][2]] |= BIT0;
                            block[index + byte_sums[i][3]] |= BIT6;
                            block[index + byte_sums[i][4]] |= BIT1;
                            block[index + byte_sums[i][5]] |= BIT7;
                            block[index + byte_sums[i][6]] |= BIT3;
                            block[index + byte_sums[i][7]] |= BIT5;
                            index += primes[i];
                        }
                        break;
                    // p == 13 mod 30
                    case 3:
                        while (index < bound) {
                            block[index]                   |= BIT3;
                            block[index + byte_sums[i][1]] |= BIT0;
                            block[index + byte_sums[i][2]] |= BIT6;
                            block[index + byte_sums[i][3]] |= BIT5;
                            block[index + byte_sums[i][4]] |= BIT2;
                            block[index + byte_sums[i][5]] |= BIT1;
                            block[index + byte_sums[i][6]] |= BIT7;
                            block[index + byte_sums[i][7]] |= BIT4;
                            index += primes[i];
                        }
                        break;
                    // p == 17 mod 30
                    case 4:
                        while (index < bound) {
                            block[index]                   |= BIT4;
                            block[index + byte_sums[i][1]] |= BIT7;
                            block[index + byte_sums[i][2]] |= BIT1;
                            block[index + byte_sums[i][3]] |= BIT2;
                            block[index + byte_sums[i][4]] |= BIT5;
                            block[index + byte_sums[i][5]] |= BIT6;
                            block[index + byte_sums[i][6]] |= BIT0;
                            block[index + byte_sums[i][7]] |= BIT3;
                            index += primes[i];
                        }
                        break;
                    // p == 19 mod 30
                    case 5:
                        while (index < bound) {
                            block[index]                   |= BIT5;
                            block[index + byte_sums[i][1]] |= BIT3;
                            block[index + byte_sums[i][2]] |= BIT7;
                            block[index + byte_sums[i][3]] |= BIT1;
                            block[index + byte_sums[i][4]] |= BIT6;
                            block[index + byte_sums[i][5]] |= BIT0;
                            block[index + byte_sums[i][6]] |= BIT4;
                            block[index + byte_sums[i][7]] |= BIT2;
                            index += primes[i];
                        }
                        break;
                    // p == 23 mod 30
                    case 6:
                        while (index < bound) {
                            block[index]                   |= BIT6;
                            block[index + byte_sums[i][1]] |= BIT2;
                            block[index + byte_sums[i][2]] |= BIT3;
                            block[index + byte_sums[i][3]] |= BIT7;
                            block[index + byte_sums[i][4]] |= BIT0;
                            block[index + byte_sums[i][5]] |= BIT4;
                            block[index + byte_sums[i][6]] |= BIT5;
                            block[index + byte_sums[i][7]] |= BIT1;
                            index += primes[i];
                        }
                        break;
                    // p == 29 mod 30
                    case 7:
                        while (index < bound) {
                            block[index]                   |= BIT7;
                            block[index + byte_sums[i][1]] |= BIT6;
                            block[index + byte_sums[i][2]] |= BIT5;
                            block[index + byte_sums[i][3]] |= BIT4;
                            block[index + byte_sums[i][4]] |= BIT3;
                            block[index + byte_sums[i][5]] |= BIT2;
                            block[index + byte_sums[i][6]] |= BIT1;
                            block[index + byte_sums[i][7]] |= BIT0;
                            index += primes[i];
                        }
                        break;
                }
                // important: put k back in it's place!
                k = 0;
            }
            while (byte_sums[i][k] < block_size - index) {
                block[index + byte_sums[i][k]] |= pattern[wheel_start[i] ][k];
                ++k;
            }
            next_byte[i] = index + offset;
            wheel_index[i] = k;
        }
        // Add up the number of zero bits.
        prime_count += count_primes(block);//, unset_bits);
        offset += block_size;
    }

    // Substract the bits that were overcounted.
    int64_t remaining_bytes = final_index - offset + block_size;
    for (int64_t i = block_size - 1; i > remaining_bytes; --i) {
        prime_count -= unset_bits[block[i]];
    }
    // Process the last byte.
    uint8_t final_byte = block[remaining_bytes];
    int64_t initial_error = unset_bits[final_byte];
    for (int64_t i = 0; i < (candidates & 7); ++i) {
        if (!(final_byte & 1)) {
            --initial_error;
        }
        final_byte >>= 1;
    }
    prime_count -= initial_error;

    // 1 was counted, but the primes 2, 3 and 5 were not.
    cout << prime_count + 7 << '\n';
    free(pre_sieve);

    return 0;
}
