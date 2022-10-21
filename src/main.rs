use std::thread;
//use std::time::Duration;
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

fn count_primes(primes: &[(usize, usize)], flags: &mut [u8]) -> u32 {
    flags.fill(0);
    sieve(primes, flags);
    flags.iter().map(|&b| b.count_zeros()).sum()
}

fn main() {
    const N: usize = 2;
    let primes = simple_sieve(96);
    dbg!(&primes);
    let mut flags: [u8; N] = [0; N];
    let x = count_primes(&[(5, 4), (7, 1)], &mut flags);
    let y = count_primes(&[(5, 0), (7, 5)], &mut flags);
    dbg!(x + y + 4);

    //let answer: u32 = count_primes(&[(2, 0), (3, 0), (5, 0), (7, 0)], &mut flags);
    /*
    for &p in &[2, 3, 5, 7] {
        sieve(p, 0, &mut flags);
    }
    for &f in &flags {
        println!("{:#b}", f);
    }
    let answer: u32 = (0..N).map(|i| flags[i].count_zeros()).sum::<u32>() + 3;
    */
    /*
    let answer: usize = (0..2).map(|i| {
        thread::spawn(move || -> usize {
            (100 * i.. 100 * (i + 1)).sum()
            //println!("hi number {total} from thread {i}");
        }).join().unwrap()
    }).sum();
    dbg!(answer);
    */
}

