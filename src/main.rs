mod sieve;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.len() {
        1 => {
            let now = std::time::Instant::now();
            let prime_count: u64 = sieve::count_primes(1);
            println!("The number of primes is {prime_count}");
            println!(
                "Time taken was {} seconds",
                now.elapsed().as_micros() as f64 / 1_000_000.
            );
        }
        2 => {
            let now = std::time::Instant::now();
            let nblocks: usize = args[1]
                .replace('_', "")
                .parse()
                .expect("Invalid numeric literal format");
            let prime_count: u64 = sieve::count_primes(nblocks);
            println!("The number of primes is {prime_count}");
            println!(
                "Time taken was {} seconds",
                now.elapsed().as_micros() as f64 / 1_000_000.
            );
        }
        _ => {
            println!("An incorrect number of arguments was given.");
        }
    }
}
