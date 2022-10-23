This is a multithreaded prime sieve written in Rust. The algorithm and it's explanation is found in src/sieve/mod.rs.
The sieve is so fast that it doesn't actually generate the primes, but only counts them, thus effectively it implements
a way to compute the prime counting function https://en.wikipedia.org/wiki/Prime-counting_function.
Although there are surely more effictive ways to obtain values of this functions much faster,
(pecially if you are only interested in an approximate answer), I believe that this prime sieve in Rust is one of
the fastest ways of computing it by 'naive' counting.

For instance, I was able to compute $\pi(27525120000) = 1197240459$ in 0.985 second on a 12th generation Core i7 cpu.
