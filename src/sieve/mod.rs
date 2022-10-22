// Takes approximately half a second to count primes below 10 billion.
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
    (11..=n)
        .filter(|&p| is_prime[p])
        .collect()
}

/// Set the flags that are k mod p in block. Return next offset.
fn sieve_with_prime(p: usize, k: usize, block: &mut [u64]) -> usize {
    let mut j = k;
    while (j >> 6) < block.len() {
        block[j >> 6] |= 1 << (j & 0x3F);
        j += p;
    }
    j - 64 * block.len()
}

/// Sieve primes from a block.
fn sieve_block(small_primes: &[usize], offsets: &mut [usize], block: &mut [u64]) {
    block.fill(0);
    small_primes
        .iter()
        .zip(offsets.iter_mut())
        .for_each(|(p, k)| *k = sieve_with_prime(*p, *k, block));
}

/// Count the zero bits in a block.
fn count_zeros(block: &[u64]) -> u32 {
    block.iter().map(|&b| b.count_zeros()).sum()
}

fn sieve_blocks(
    nblocks: usize,
    small_primes: &[usize],
    offsets: &mut [usize],
    block: &mut [u64],
) -> u64 {
    (0..nblocks)
        .map(|_| {
            sieve_block(small_primes, offsets, block);
            count_zeros(block) as u64
        })
        .sum()
}

/// Compute a / b mod WHEEL_SZ
fn mod_div(a: usize, b: usize) -> Option<usize> {
    WHEEL.iter().find(|&&w| (w * b) % WHEEL_SZ == a).copied()
}

/// Count the number of primes congruent to w modulo WHEEL_SZ.
fn count_res_class(small_primes: &[usize], nblocks: usize, w: usize) -> u64 {
    let mut block: [u64; BLOCK_SZ] = [0; BLOCK_SZ];
    let mut offsets = small_primes
        .iter()
        .map(|&p| (mod_div(w, p).unwrap() * p) / WHEEL_SZ)
        .collect::<Vec<usize>>();
    sieve_blocks(nblocks, small_primes, &mut offsets, &mut block)
}

pub fn count_primes(nblocks: usize) -> u64 {
    let n: usize = nblocks * BLOCK_SZ * WHEEL_SZ * 64;
    println!("Computing primes below {n}");
    let s = (n as f64).sqrt() as usize;
    let small_primes: Vec<usize> = small_sieve(s);

    let handles = WHEEL
        .iter()
        .map(|&w| {
            let small_primes = small_primes.clone();
            std::thread::spawn(move || -> u64 { count_res_class(&small_primes, nblocks, w) })
        })
        .collect::<Vec<_>>();

    let answer: u64 = handles
        .into_iter()
        .map(|handle| handle.join().unwrap())
        .sum();

    // Add the primes in the wheel that were not tested.
    // Also add primes used for sieving because they sieved themselves out.
    // Finally subtract 1 because it isn't prime.
    answer + small_primes.len() as u64 + 3
}
