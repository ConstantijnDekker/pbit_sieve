/// Algorithm description:
/// For each residue w in WHEEL, we start a thread to mark off composite numbers
/// that are congruent to w mod WHEEL_SZ. We do this by dividing these
/// candidates in groups of BLOCK_SZ WORD_SZ-bit words, for cache efficiency,
/// with the lth bit of the mth-block representing
///     (BLOCK_SZ * WORD_SZ * m + l) * WHEEL_SZ + w.
///
/// We know that for a given prime p not in WHEEL_PRIMES there exists a unique
/// r in WHEEL such that p * r == w mod WHEEL_SZ,
/// and all solutions to this equation are given by the set
///     {r + WHEEL_SZ * l : l an integer}.
/// Because we have the equality
///     (p * (r + WHEEL_SZ * l)) / WHEEL_SZ == (p * r) / WHEEL_SZ + l * p
/// we can mark off such multiples by moving p bits ahead.

const WORD_SZ: usize = 64;
const WHEEL_PRIMES: [usize; 4] = [2, 3, 5, 7];
const WHEEL: [usize; 48] = [
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181,
    187, 191, 193, 197, 199, 209,
];
const WHEEL_SZ: usize = 210;
const BLOCK_SZ: usize = 1 << 12;

/// Plain sieve of Eratosthenes that finds all primes <= n that are not in the WHEEL.
fn small_sieve(n: usize) -> Vec<usize> {
    let mut is_prime = vec![true; n + 1];

    for i in (2..=n).take_while(|i| i * i <= n) {
        if is_prime[i] {
            for j in (i * i..=n).step_by(i) {
                is_prime[j] = false;
            }
        }
    }

    (2..=n)
        .filter(|&p| is_prime[p] && !WHEEL_PRIMES.contains(&p))
        .collect()
}

/// Set the flags that are k plus a multiple of p in block. Return next offset.
fn sieve_with_prime(p: usize, k: usize, block: &mut [u64]) -> usize {
    let mut j = k;
    while j < WORD_SZ * block.len() {
        block[j / WORD_SZ] |= 1 << (j % WORD_SZ);
        j += p;
    }
    j - WORD_SZ * block.len()
}

/// Sieve primes from a block.
fn sieve_block(primes: &[usize], offsets: &mut [usize], block: &mut [u64]) {
    block.fill(0);
    primes
        .iter()
        .zip(offsets.iter_mut())
        .for_each(|(p, k)| *k = sieve_with_prime(*p, *k, block));
}

/// Count the zero bits in a block.
fn count_zeros(block: &[u64]) -> u64 {
    block.iter().map(|&b| b.count_zeros() as u64).sum()
}

/// Return bitmask for lowest 'width' bits of a word.
fn mask(width: usize) -> u64 {
    (1u64 << width).wrapping_sub(1)
}

/// Sieve nbits amount of bits with the primes at the given offsets.
fn sieve(mut nbits: usize, primes: &[usize], offsets: &mut [usize], block: &mut [u64]) -> u64 {
    let mut prime_count = 0;
    while nbits >= WORD_SZ * BLOCK_SZ {
        sieve_block(primes, offsets, block);
        prime_count += count_zeros(block);
        nbits -= WORD_SZ * BLOCK_SZ;
    }
    // The final block does not have to be fully processed.
    let last = nbits / WORD_SZ;
    sieve_block(primes, offsets, &mut block[..=last]);
    prime_count += count_zeros(&mut block[..=last]);

    // small correction for overcounting zeros in the last word
    prime_count - (block[last] | mask(nbits % WORD_SZ)).count_zeros() as u64
}

/// Compute a / b mod WHEEL_SZ
fn mod_div(a: usize, b: usize) -> Option<usize> {
    WHEEL.iter().find(|&&w| (w * b) % WHEEL_SZ == a).copied()
}

/// Count the number of primes congruent to w modulo WHEEL_SZ.
fn count_res_class(nbits: usize, primes: &[usize], w: usize) -> u64 {
    let mut block: [u64; BLOCK_SZ] = [0; BLOCK_SZ];
    // Offsets[i] is the bitindex primes[i] has to start
    // marking off in the next block.
    let mut offsets = primes
        .iter()
        .map(|&p| -> usize {
            let mut r = mod_div(w, p).unwrap();
            // Make sure r >= p while remaining the same mod WHEEL_SZ
            if r >= (p % WHEEL_SZ) {
                r += p - (p % WHEEL_SZ)
            } else {
                r += p - (p % WHEEL_SZ) + WHEEL_SZ
            };
            // bit index of first multiple of p to be crossed of
            (r * p) / WHEEL_SZ
        })
        .collect::<Vec<usize>>();
    sieve(nbits, primes, &mut offsets, &mut block)
}

/// Compute total number of primes <= n.
pub fn count_primes(n: usize) -> u64 {
    // Sieve with the primes below the square root of n.
    let primes: Vec<usize> = small_sieve((n as f64).sqrt() as usize);

    // Start a thread for each residue in our wheel
    let handles = WHEEL
        .iter()
        .map(|&w| {
            let primes = primes.clone();
            std::thread::spawn(move || -> u64 {
                if w <= (n % WHEEL_SZ) {
                    count_res_class(n / WHEEL_SZ + 1, &primes, w)
                } else {
                    count_res_class(n / WHEEL_SZ, &primes, w)
                }
            })
        })
        .collect::<Vec<_>>();

    let mut prime_count: u64 = handles
        .into_iter()
        .map(|handle| handle.join().unwrap())
        .sum();

    // Add the primes used for the WHEEL.
    prime_count += WHEEL_PRIMES.iter().filter(|&&p| p <= n).count() as u64;

    // Finally subtract 1 because it isn't prime.
    if n >= 1 {
        prime_count -= 1;
    }

    prime_count
}
