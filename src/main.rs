mod sieve;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    match args.len() {
        2 => {
            let now = std::time::Instant::now();
            let n: usize = args[1]
                .replace('_', "")
                .parse()
                .expect("Invalid numeric literal format");
            let prime_count: u64 = sieve::count_primes(n);
            println!("Found {prime_count} primes below {n} in {} seconds",
                now.elapsed().as_micros() as f64 / 1_000_000.
            );
        }
        _ => {
            println!("Usage: PROGRAM_NAME [limit]");
        }
    }
}
