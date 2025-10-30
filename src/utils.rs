use bio_seq::seq::Seq;
use statrs::distribution::{Binomial, DiscreteCDF};
use memory_stats::memory_stats;

pub fn log_memory_usage(info: bool, message: &str) {
    if let Some(usage) = memory_stats() {
        if info{
            log::info!(
                "{} --- Memory usage: {:.2} GB",
                message,
                usage.physical_mem as f64 / 1_000_000_000.
            );
        }
        else{
            log::debug!(
                "{} --- Memory usage: {:.2} GB",
                message,
                usage.physical_mem as f64 / 1_000_000_000.
            );
        }
    }
    else{
        log::info!("Memory usage: unknown (WARNING)");
    }
}

#[inline]
pub fn div_rounded(a: usize, b: usize) -> usize {
    (a + b / 2) / b
}



pub fn first_word(s: &str) -> String{
    s.split_whitespace().next().unwrap_or(s).to_string()
}

pub fn binomial_test(n: u64, k: u64, p: f64) -> f64 {
    // n: number of trials
    // k: number of successes
    // p: probability of success

    // Create a binomial distribution
    let binomial = Binomial::new(p, n).unwrap();

    // Calculate the probability of observing k or more successes
    let p_value = 1.0 - binomial.cdf(k);

    p_value
}

pub fn reverse_complement(seq: &[u8]) -> Vec<u8> {
    let mut revcomp = Vec::with_capacity(seq.len());
    for &base in seq.iter().rev() {
        let comp_base = match base {
            b'A' | b'a' => b'T',
            b'T' | b't' => b'A',
            b'C' | b'c' => b'G',
            b'G' | b'g' => b'C',
            b'N' | b'n' => b'N',
            _ => b'N', // Handle unexpected characters
        };
        revcomp.push(comp_base);
    }
    revcomp
}
