const WHEEL_PRIMES: [usize; 4] = [2, 3, 5, 7];
const WHEEL: [usize; 48] = [
    1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101,
    103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181,
    187, 191, 193, 197, 199, 209,
];
const WHEEL_SZ: usize = 210;
const BLOCK_SZ: usize = 1 << 15;

/// Plain sieve of Eratosthenes that finds all primes <= n that are not in the WHEEL.
fn simple_sieve(n: usize) -> Vec<usize> {
    let mut is_prime = vec![true; n + 1];
    for i in (2..=n).take_while(|i| i * i <= n) {
        for j in (i * i..=n).step_by(i) {
            is_prime[j] = false;
        }
    }
    (2..=n)
        .filter(|&p| is_prime[p] && !WHEEL_PRIMES.contains(&p))
        .collect()
}

/// Set the flags that are k mod p in block. Return next offset.
fn sieve_with_prime(p: usize, k: usize, block: &mut [u8]) -> usize {
    let mut j = k;
    while (j >> 3) < block.len() {
        block[j >> 3] |= 1 << (j & 0b111);
        j += p;
    }
    j - 8 * block.len()
}

/// Sieve primes from a block.
fn sieve_block(records: &mut [(usize, usize)], block: &mut [u8]) {
    block.fill(0);
    records
        .iter_mut()
        .for_each(|(p, k)| *k = sieve_with_prime(*p, *k, block));
}

/// Count the zero bits in a block.
fn count_zeros(block: &[u8]) -> u32 {
    block.iter().map(|&b| b.count_zeros()).sum()
}

fn sieve_blocks(nblocks: usize, records: &mut [(usize, usize)], block: &mut [u8]) -> u32 {
    (0..nblocks)
        .map(|_| {
            sieve_block(records, block);
            count_zeros(block)
        })
        .sum()
}

/// Compute a / b mod WHEEL_SZ
fn mod_div(a: usize, b: usize) -> Option<usize> {
    WHEEL.iter().find(|&&w| (w * b) % WHEEL_SZ == a).copied()
}

/// Count the number of primes congruent to w modulo WHEEL_SZ.
fn count_res_class(sieving_primes: Vec<usize>, nblocks: usize, w: usize) -> u32 {
    let mut block: [u8; BLOCK_SZ] = [0; BLOCK_SZ];
    let mut records = sieving_primes
        .iter()
        .map(|&p| (p, (mod_div(w, p).unwrap() * p) / WHEEL_SZ))
        .collect::<Vec<(usize, usize)>>();
    sieve_blocks(nblocks, &mut records, &mut block)
}

pub fn count_primes(nblocks: usize) -> u64 {
    let byte_sz = 8 * WHEEL_SZ;
    let n: usize = nblocks * BLOCK_SZ * byte_sz;
    println!("Computing primes below {n}");
    let s = (n as f64).sqrt() as usize;
    let sieving_primes: Vec<usize> = simple_sieve(s);

    let handles = WHEEL
        .iter()
        .map(|&w| {
            let sieving_primes = sieving_primes.clone();
            std::thread::spawn(move || -> u32 { count_res_class(sieving_primes, nblocks, w) })
        })
        .collect::<Vec<_>>();

    let answer: u32 = handles
        .into_iter()
        .map(|handle| handle.join().unwrap())
        .sum();

    // Add the primes in the wheel that were not tested.
    // Also add primes used for sieving because they sieved themselves out.
    // Finally subtract 1 because it isn't prime.
    answer as u64 + WHEEL_PRIMES.len() as u64 + sieving_primes.len() as u64 - 1
}
