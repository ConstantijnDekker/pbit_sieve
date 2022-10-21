/*
const WHEEL_PRIMES: [usize; 2] = [2, 3];
const WHEEL: [usize; 2] = [1, 5];
const WHEEL_SZ: usize = 6;
*/


//const WHEEL_PRIMES: [usize; 3] = [2, 3, 5];
//const WHEEL: [usize; 8] = [1, 7, 11, 13, 17, 19, 23, 29];
//const WHEEL_SZ: usize = 30;

const WHEEL_PRIMES: [usize; 4] = [2, 3, 5, 7];
const WHEEL: [usize; 48] = [1, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47, 53, 59, 61, 67, 71, 73, 79, 83, 89, 97, 101, 103, 107, 109, 113, 121, 127, 131, 137, 139, 143, 149, 151, 157, 163, 167, 169, 173, 179, 181, 187, 191, 193, 197, 199, 209];
const WHEEL_SZ: usize = 210;

const BLOCK_SZ: usize = 1 << 16;

fn simple_sieve(n: usize) -> Vec<usize> {
    let mut is_prime = vec![true; n + 1];
    for i in (2..=n).take_while(|i| i * i <= n) {
        for j in (i * i..=n).step_by(i) {
            is_prime[j] = false;
        }
    }
    (2..=n).filter(|&p| is_prime[p] && !WHEEL_PRIMES.contains(&p)).collect()
}

fn set_bit(flags: &mut [u8], bit_idx: usize) {
    flags[bit_idx >> 3] |= 1 << (bit_idx & 0b111);
}

// Set flags that are k mod p in flags. Return next offset.
fn sieve_with_prime(p: usize, k: usize, flags: &mut [u8]) -> usize {
    let mut j = k;
    while j < 8 * flags.len() {
        set_bit(flags, j);
        j += p;
    }
    j - 8 * flags.len()
}

fn sieve(records: &mut [(usize, usize)], flags: &mut [u8]) {
    flags.fill(0);
    records.iter_mut().for_each(|(p, k)|
        *k = sieve_with_prime(*p, *k, flags));
}

fn count_zeros(flags: &[u8]) -> u32 {
    flags.iter().map(|&b| b.count_zeros()).sum()
}

fn sieve_blocks(num_iter: usize, records: &mut [(usize, usize)], flags: &mut [u8]) -> u32 {
    (0..num_iter).map(|_| {
            sieve(records, flags);
            count_zeros(flags)
    }).sum()
}

// Compute a / b mod WHEEL_SZ
fn mod_div(a: usize, b: usize) -> Option<usize> {
    WHEEL.iter().find(|&&w| (w * b) % WHEEL_SZ == a).copied()
}

fn count_primes(small_primes: Vec<usize>, nblocks: usize, w: usize) -> u32 {
    let mut flags = vec![0; BLOCK_SZ];
    let mut records = small_primes.iter().map(|&p| {
        (p, (mod_div(w, p).unwrap() * p) / WHEEL_SZ)
    }).collect::<Vec<(usize, usize)>>();
    sieve_blocks(nblocks, &mut records, &mut flags)
}

fn main() {
    let byte_sz = 8 * WHEEL_SZ;
    let nblocks = 9;
    let n: usize = nblocks * BLOCK_SZ * byte_sz;
    println!("Calculating number of primes below {n}");
    let s = (n as f64).sqrt() as usize;
    let small_primes: Vec<usize> = simple_sieve(s);
    let now = std::time::Instant::now();

    /*
    // Single threaded
    let answer: u32 = WHEEL.iter()
        .map(|&w| count_primes(small_primes.clone(), nbytes, w))
        .sum();
    */
    let handles = WHEEL.iter().map(|&w| {
        let small_primes = small_primes.clone();
        std::thread::spawn(move || -> u32 {
            count_primes(small_primes, nblocks, w)
        })
    }).collect::<Vec<_>>();
    let answer: u32 = handles.into_iter().map(|handle| handle.join().unwrap()).sum();
    // Add the primes in the wheel that were not tested.
    // Also add primes used for sieving because they sieved themselves out.
    // Finally subtract 1 because it isn't prime.
    let final_answer = answer as usize + WHEEL_PRIMES.len() + small_primes.len() - 1;
    dbg!(final_answer);
    println!("Time taken was {} milliseconds", now.elapsed().as_micros());
}


