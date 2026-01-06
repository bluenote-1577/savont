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

/// Homopolymer compress a sequence
/// Returns (hpc_sequence, hp_lengths) where hp_lengths[i] is the run length of hpc_sequence[i]
/// Example: b"AAACGT" -> (b"ACGT", vec![3, 1, 1, 1])
pub fn homopolymer_compress(seq: &[u8], do_hpc: bool) -> (Vec<u8>, Vec<u8>) {
    if seq.is_empty() {
        return (Vec::new(), Vec::new());
    }

    let mut hpc_seq = Vec::with_capacity(seq.len());
    let mut hp_lengths = Vec::with_capacity(seq.len());

    let mut current_base = seq[0];
    let mut current_length = 1u8;
    if do_hpc{
        for i in 1..seq.len() {
            if seq[i] == current_base && current_length < 255 {
                // Continue the homopolymer run (cap at 255)
                current_length += 1;
            } else {
                // End of run, save it
                hpc_seq.push(current_base);
                hp_lengths.push(current_length);

                // Start new run
                current_base = seq[i];
                current_length = 1;
            }
        }
    }
    else{
        hpc_seq.extend_from_slice(seq);
        hp_lengths.extend(vec![1u8; seq.len()]);
    }

    // Don't forget the last run
    hpc_seq.push(current_base);
    hp_lengths.push(current_length);

    hpc_seq.shrink_to_fit();
    hp_lengths.shrink_to_fit();

    (hpc_seq, hp_lengths)
}

/// Decompress a homopolymer-compressed sequence
/// Takes (hpc_sequence, hp_lengths) and returns the full sequence
/// Example: (b"ACGT", vec![3, 1, 1, 1]) -> b"AAACGT"
pub fn homopolymer_decompress(hpc_seq: &[u8], hp_lengths: &[u8]) -> Vec<u8> {
    if hpc_seq.len() != hp_lengths.len() {
        log::warn!("HPC sequence and lengths mismatch: {} vs {}", hpc_seq.len(), hp_lengths.len());
        return hpc_seq.to_vec(); // Return as-is if mismatch
    }

    let total_length: usize = hp_lengths.iter().map(|&l| l as usize).sum();
    let mut seq = Vec::with_capacity(total_length);

    for (&base, &length) in hpc_seq.iter().zip(hp_lengths.iter()) {
        for _ in 0..length {
            seq.push(base);
        }
    }

    seq
}

/// Homopolymer compress a sequence with quality scores
/// Returns (hpc_sequence, hpc_qualities, hp_lengths)
/// Quality strategy: use minimum quality from each homopolymer run (most conservative)
/// Example: (b"AAACGT", [30,35,40,25,30,35]) -> (b"ACGT", [30,25,30,35], [3,1,1,1])
pub fn homopolymer_compress_with_quality(seq: &[u8], qualities: &[u8], do_hpc: bool) -> (Vec<u8>, Vec<u8>, Vec<u8>) {
    if seq.is_empty() || seq.len() != qualities.len() {
        return (Vec::new(), Vec::new(), Vec::new());
    }

    let mut hpc_seq = Vec::with_capacity(seq.len());
    let mut hpc_qualities = Vec::with_capacity(seq.len());
    let mut hp_lengths = Vec::with_capacity(seq.len());

    let mut current_base = seq[0];
    let mut current_length = 1u8;
    let mut current_min_quality = qualities[0];

    if do_hpc{
        for i in 1..seq.len() {
            if seq[i] == current_base && current_length < 255 {
                // Continue the homopolymer run
                current_length += 1;
                current_min_quality = current_min_quality.min(qualities[i]);
            } else {
                // End of run, save it
                hpc_seq.push(current_base);
                hpc_qualities.push(current_min_quality);
                hp_lengths.push(current_length);

                // Start new run
                current_base = seq[i];
                current_length = 1;
                current_min_quality = qualities[i];
            }
        }

        // Don't forget the last run
        hpc_seq.push(current_base);
        hpc_qualities.push(current_min_quality);
        hp_lengths.push(current_length);
    }
    else{
        hpc_seq.extend_from_slice(seq);
        hpc_qualities.extend_from_slice(qualities);
        hp_lengths.extend(vec![1u8; seq.len()]);
    }

    hpc_seq.shrink_to_fit();
    hpc_qualities.shrink_to_fit();
    hp_lengths.shrink_to_fit();

    (hpc_seq, hpc_qualities, hp_lengths)
}

/// Expand binned quality scores to match sequence length
/// Takes a quality sequence iterator, multiplies by 3, adds 33 offset,
/// expands by bin_size, and adjusts to match sequence length
pub fn expand_binned_qualities_from_iter<I>(qual_iter: I, seq_len: usize, bin_size: usize) -> Vec<u8>
where
    I: Iterator<Item = u8>,
{
    // Multiply quality by 3 and add 33 (Phred+33 offset)
    let qual_u8: Vec<u8> = qual_iter.map(|x| (x as u8) * 3 + 33).collect();

    // Expand each quality by bin_size
    let mut expanded_quals: Vec<u8> = qual_u8
        .iter()
        .flat_map(|&x| vec![x; bin_size])
        .collect();

    // Adjust to match sequence length
    if expanded_quals.len() > seq_len {
        expanded_quals.truncate(seq_len);
    } else if expanded_quals.len() < seq_len {
        let last_qual = *expanded_quals.last().unwrap_or(&33);
        expanded_quals.extend(vec![last_qual; seq_len - expanded_quals.len()]);
    }

    expanded_quals
}
