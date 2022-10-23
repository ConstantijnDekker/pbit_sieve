const WHEEL: [usize; 48] = [
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181,
    187, 191, 193, 197, 199, 209,
];
const WHEEL_SZ: usize = 210;
const BLOCK_SZ: usize = 1 << 12;

/// Plain sieve of Eratosthenes that finds all primes <= n that are
/// not in the WHEEL.
fn small_sieve(n: usize) -> Vec<usize> {
    let mut is_prime = vec![true; n + 1];

    for i in (2..=n).take_while(|i| i * i <= n) {
        if is_prime[i] {
            for j in (i * i..=n).step_by(i) {
                is_prime[j] = false;
            }
        }
    }
    // Start from 11 to exclude primes 2, 3, 5, 7
    (11..=n).filter(|&p| is_prime[p]).collect()
}

/// Set the flags that are k mod p in block. Return next offset.
fn sieve_with_prime(p: usize, k: usize, block: &mut [u64]) -> usize {
    let mut j = k;
    while j < 64 * block.len() {
        block[j >> 6] |= 1 << (j & 0x3F);
        j += p;
    }
    j - 64 * block.len()
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
    while nbits >= 64 * BLOCK_SZ {
        sieve_block(primes, offsets, block);
        prime_count += count_zeros(block);
        nbits -= 64 * BLOCK_SZ;
    }
    // The final block does not have to be fully processed.
    let last = nbits / 64;
    sieve_block(primes, offsets, &mut block[..=last]);
    prime_count += count_zeros(&mut block[..=last]);

    // small correction for overcounting zeros in the last word
    prime_count - (block[last] | mask(nbits % 64)).count_zeros() as u64
}

/// Compute a / b mod WHEEL_SZ
fn mod_div(a: usize, b: usize) -> Option<usize> {
    WHEEL.iter().find(|&&w| (w * b) % WHEEL_SZ == a).copied()
}

/// Count the number of primes congruent to w modulo WHEEL_SZ.
fn count_res_class(nbits: usize, primes: &[usize], w: usize) -> u64 {
    let mut block: [u64; BLOCK_SZ] = [0; BLOCK_SZ];
    let mut offsets = primes
        .iter()
        .map(|&p| -> usize {
            let mut r = mod_div(w, p).unwrap();
            if r >= (p % WHEEL_SZ) {
                r += p - (p % WHEEL_SZ)
            } else {
                r += p - (p % WHEEL_SZ) + WHEEL_SZ
            };
            (r * p) / WHEEL_SZ
        })
        .collect::<Vec<usize>>();
    sieve(nbits, primes, &mut offsets, &mut block)
}

/// Compute total number of primes <= n.
pub fn count_primes(n: usize) -> u64 {
    let s = (n as f64).sqrt() as usize;
    let primes: Vec<usize> = small_sieve(s);
    let nbits = n / WHEEL_SZ;

    // Start a thread for each residue in our wheel
    let handles = WHEEL
        .iter()
        .map(|&w| {
            let primes = primes.clone();
            std::thread::spawn(move || -> u64 {
                if w <= (n % WHEEL_SZ) {
                    count_res_class(nbits + 1, &primes, w)
                } else {
                    count_res_class(nbits, &primes, w)
                }
            })
        })
        .collect::<Vec<_>>();

    let mut prime_count: u64 = handles
        .into_iter()
        .map(|handle| handle.join().unwrap())
        .sum();

    // Add the primes in the wheel that were not tested.
    for &p in &[2, 3, 5, 7] {
        if n >= p {
            prime_count += 1;
        }
    }
    // Finally subtract 1 because it isn't prime.
    if n >= 1 {
        prime_count -= 1;
    }

    prime_count
}
