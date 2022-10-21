use std::thread;
//use std::time::Duration;
const WHEEL_PRIMES: [usize; 2] = [2, 3];
const WHEEL: [usize; 2] = [1, 5];
const WHEEL_SZ: usize = 6;

fn simple_sieve(n: usize) -> Vec<usize> {
    let mut is_prime = vec![true; n + 1];
    for i in (2..=n).take_while(|i| i * i <= n) {
        for j in (i * i..=n).step_by(i) {
            is_prime[j] = false;
        }
    }
    (2..=n).filter(|&p| is_prime[p]).collect()
}

fn set_bit(flags: &mut [u8], bit_idx: usize) {
    flags[bit_idx >> 3] |= 1 << (bit_idx & 0b111);
}

fn sieve_with_prime(p: usize, k: usize, flags: &mut [u8]) {
    for j in (k..flags.len() * 8).step_by(p) {
        set_bit(flags, j);
    }
}

fn sieve(primes: &[(usize, usize)], flags: &mut [u8]) {
    primes.iter().for_each(|&(p, k)| sieve_with_prime(p, k, flags));
}

fn count_bits(primes: &[(usize, usize)], flags: &mut [u8]) -> u32 {
    flags.fill(0);
    sieve(primes, flags);
    flags.iter().map(|&b| b.count_zeros()).sum()
}

fn mod_div(a: usize, b: usize) -> Option<usize> {
    WHEEL.iter().find(|&&w| (w * b) % WHEEL_SZ == a).copied()
}

fn count_primes(small_primes: Vec<usize>, nbytes: usize, w: usize) -> u32 {
    let mut flags = vec![0; nbytes];
    let records = small_primes.iter().map(|&p| {
        (p, (mod_div(w, p).unwrap() * p) / WHEEL_SZ)
    }).collect::<Vec<(usize, usize)>>();
    count_bits(&records, &mut flags)
}

fn main() {
    let byte_sz = 8 * WHEEL_SZ;
    let n: usize = 1_000_000 * byte_sz;
    let nbytes = n / byte_sz;
    let s = (n as f64).sqrt() as usize;
    let small_primes: Vec<usize> = simple_sieve(s);
    let small_primes: Vec<usize> = small_primes.iter()
        .filter(|p| !WHEEL_PRIMES.contains(p))
        .cloned()
        .collect::<Vec<_>>();
    let now = std::time::Instant::now();

    /*
    // Single threaded
    let answer: u32 = WHEEL.iter()
        .map(|&w| count_primes(small_primes.clone(), nbytes, w))
        .sum();
    // Multithreaded 1: slow
    let answer: u32 = WHEEL.iter().map(|&w| {
        let small_primes = small_primes.clone();
        thread::spawn(move || -> u32 {
            count_primes(small_primes, nbytes, w)
        }).join().unwrap()
    }).sum();
    */
    // Multithreaded 2: fast
    let handles = WHEEL.iter().map(|&w| {
        let small_primes = small_primes.clone();
        thread::spawn(move || -> u32 {
            count_primes(small_primes, nbytes, w)
        })
    }).collect::<Vec<_>>();
    let answer: u32 = handles.into_iter().map(|handle| handle.join().unwrap()).sum();
    let final_answer = answer as usize + WHEEL_PRIMES.len() + small_primes.len() - 1;
    dbg!(final_answer);
    println!("Time taken was {} milliseconds", now.elapsed().as_micros());
}


