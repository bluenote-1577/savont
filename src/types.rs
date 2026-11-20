//Various byte-tables and hashing methods are taken from miniprot by Heng Li. Attached below is their license:
//The MIT License

// **** miniprot LICENSE ***
//Copyright (c) 2022-     Dana-Farber Cancer Institute
//
//Permission is hereby granted, free of charge, to any person obtaining
//a copy of this software and associated documentation files (the
//"Software"), to deal in the Software without restriction, including
//without limitation the rights to use, copy, modify, merge, publish,
//distribute, sublicense, and/or sell copies of the Software, and to
//permit persons to whom the Software is furnished to do so, subject to
//the following conditions:
//
//The above copyright notice and this permission notice shall be
//included in all copies or substantial portions of the Software.
//
//THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
//EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
//MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
//NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
//BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
//ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
//CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//SOFTWARE.
//******************************

use bio_seq::prelude::*;
use fxhash::FxHashMap;
use fxhash::FxHashSet;
use serde::{Deserialize, Serialize};
use smallvec::SmallVec;
use std::cmp::Ordering;
use std::collections::BTreeMap;
use std::collections::HashMap;
use std::collections::HashSet;
use std::hash::Hash;
use std::hash::{BuildHasherDefault, Hasher};
use std::path::PathBuf;

use crate::constants::ID_THRESHOLD_ITERS;
use crate::constants::MAX_GAP_CHAINING;
use crate::constants::{LSH_BUCKET_SIZE, LSH_NUM_TABLES};

pub type NodeMap<K, V> = BTreeMap<K, V>;
pub type Kmer64 = u64;
pub type Kmer32 = u32;
pub type KmerHash64 = u64;
pub type KmerHash32 = u32;
pub type Snpmer64 = u64;
pub type Splitmer64 = u64;
pub type MultiCov = [f64; ID_THRESHOLD_ITERS];

/// Represents a blockmer with its sequence and orientation
/// Structure: [anchor k-mer (k bases)][suffix (l bases)]
/// is_forward indicates whether the anchor is at the start (true) or end (false) of the blockmer
#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize, Deserialize)]
pub struct Blockmer {
    /// The full blockmer sequence (anchor + suffix) encoded in 2 bits per base
    pub kmer: u64,
    /// True if this blockmer is from forward strand (anchor at start, suffix at end)
    /// False if from reverse strand (suffix at start after RC, anchor at end)
    pub is_forward: bool,
}

impl Blockmer {
    pub fn new(kmer: u64, is_forward: bool) -> Self {
        Self { kmer, is_forward }
    }
}

impl Ord for Blockmer {
    fn cmp(&self, other: &Self) -> std::cmp::Ordering {
        // First compare by kmer, then by orientation
        self.kmer
            .cmp(&other.kmer)
            .then(self.is_forward.cmp(&other.is_forward))
    }
}

impl PartialOrd for Blockmer {
    fn partial_cmp(&self, other: &Self) -> Option<std::cmp::Ordering> {
        Some(self.cmp(other))
    }
}

impl Hash for Blockmer {
    fn hash<H: std::hash::Hasher>(&self, state: &mut H) {
        self.kmer.hash(state);
    }
}

pub const BYTE_TO_SEQ: [u8; 256] = [
    0, 1, 2, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 1, 0, 0, 0, 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 3, 3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
    0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0,
];

#[derive(Debug, Clone, Copy, PartialEq, Eq, Hash, Serialize, Deserialize)]
pub struct Kmer48 {
    data: [u8; 6],
}

impl Kmer48 {
    // Be careful about endian
    #[inline]
    pub fn from_u64(n: u64) -> Self {
        let bytes = n.to_le_bytes();

        debug_assert!(bytes[6..8].iter().all(|&b| b == 0));

        Self {
            data: [bytes[0], bytes[1], bytes[2], bytes[3], bytes[4], bytes[5]],
        }
    }

    // Assume Kmer48 is stored as little endian. So AAA = ...0000000_11_11_11
    pub fn to_u64(self) -> u64 {
        let mut bytes = [0; 8];
        bytes[0..6].copy_from_slice(&self.data);
        u64::from_le_bytes(bytes)
    }
}

impl From<u64> for Kmer48 {
    fn from(value: u64) -> Self {
        Kmer48::from_u64(value)
    }
}

/// Consensus polymorphic marker statistics
/// Contains the consensus k-mer for a polymorphic marker (SNPmer or blockmer) along with metadata
#[derive(Debug, Clone, Copy, PartialEq, Eq)]
pub struct ConsensusPoly {
    /// Median position along the read where this polymorphic marker occurs
    pub position: u32,
    /// The splitmer (anchor k-mer): masked k-mer for SNPmers, anchor k-mer for blockmers
    pub splitmer: u64,
    /// The consensus full k-mer (most common variant): full SNPmer or full blockmer
    pub kmer: Kmer48,
    /// Number of times this consensus k-mer was observed
    pub count: u32,
}

impl ConsensusPoly {
    pub fn new(position: u32, splitmer: u64, kmer: Kmer48, count: u32) -> Self {
        Self {
            position,
            splitmer,
            kmer,
            count,
        }
    }
}

// Type alias for backwards compatibility
pub type ConsensusSnpmer = ConsensusPoly;

/// Consensus sequence with depth information
/// Contains the consensus sequence and the number of reads used to generate it
#[derive(Debug, Clone, PartialEq, Eq, Default)]
pub struct ConsensusSequence {
    /// The HPC consensus sequence
    pub sequence: Vec<u8>,
    /// The homopolymer lengths corresponding to each base in sequence
    pub hp_lengths: Vec<u8>,
    /// The decompressed sequence (for merging and chimera detection)
    pub decompressed_sequence: Option<Vec<u8>>,
    /// The depth (number of reads) used to generate this consensus
    pub depth: usize,

    pub appended_depth: usize,

    pub cluster: Vec<usize>,

    pub low_quality_positions: Vec<usize>,

    pub chimera_score: Option<i64>,

    pub id: usize,

    pub unambig_best_read_map_count: Option<usize>,

    pub ambig_read_map_count: Option<usize>,

    pub num_map_leq_10nm: Option<usize>,

    /// Per-sample depths from --pooled-samples mode; empty in single-sample mode.
    pub per_sample_depths: Vec<usize>,
}

impl ConsensusSequence {
    pub fn new(
        sequence: Vec<u8>,
        hp_lengths: Vec<u8>,
        depth: usize,
        id: usize,
        cluster: Vec<usize>,
    ) -> Self {
        Self {
            sequence,
            hp_lengths,
            decompressed_sequence: None,
            depth,
            appended_depth: 0,
            low_quality_positions: Vec::new(),
            cluster,
            chimera_score: None,
            id,
            unambig_best_read_map_count: None,
            ambig_read_map_count: None,
            num_map_leq_10nm: None,
            per_sample_depths: Vec::new(),
        }
    }

    /// Decompress the HPC sequence and store it
    pub fn decompress(&mut self) {
        self.decompressed_sequence = Some(crate::utils::homopolymer_decompress(
            &self.sequence,
            &self.hp_lengths,
        ));
        let first_no_n = self
            .decompressed_sequence
            .as_ref()
            .unwrap()
            .iter()
            .position(|&b| b != b'N')
            .unwrap_or(0);
        let last_no_n = self
            .decompressed_sequence
            .as_ref()
            .unwrap()
            .iter()
            .rposition(|&b| b != b'N')
            .unwrap_or(self.decompressed_sequence.as_ref().unwrap().len() - 1);
        self.decompressed_sequence =
            Some(self.decompressed_sequence.as_ref().unwrap()[first_no_n..=last_no_n].to_vec());
    }

    /// Get the decompressed sequence, decompressing if needed
    pub fn get_decompressed(&mut self) -> &Vec<u8> {
        if self.decompressed_sequence.is_none() {
            self.decompress();
        }
        self.decompressed_sequence.as_ref().unwrap()
    }
}

#[inline]
pub fn mm_hash_64(key: u64) -> usize {
    let mut key = key;
    key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key as usize;
}

#[inline]
pub fn mm_hash(bytes: &[u8]) -> usize {
    let mut key = usize::from_ne_bytes(bytes.try_into().unwrap()) as usize;
    key = (!key).wrapping_add(key << 21); // key = (key << 21) - key - 1;
    key = key ^ key >> 24;
    key = (key.wrapping_add(key << 3)).wrapping_add(key << 8); // key * 265
    key = key ^ key >> 14;
    key = (key.wrapping_add(key << 2)).wrapping_add(key << 4); // key * 21
    key = key ^ key >> 28;
    key = key.wrapping_add(key << 31);
    return key;
}

pub struct MMHasher {
    hash: usize,
}

impl Hasher for MMHasher {
    #[inline]
    fn write(&mut self, bytes: &[u8]) {
        self.hash = mm_hash(bytes);
    }
    #[inline]
    fn finish(&self) -> u64 {
        self.hash as u64
    }
}

impl Default for MMHasher {
    #[inline]
    fn default() -> MMHasher {
        MMHasher { hash: 0 }
    }
}

//Implement minimap2 hashing, will test later.
pub type MMBuildHasher = BuildHasherDefault<MMHasher>;
pub type MMHashMap<K, V> = HashMap<K, V, MMBuildHasher>;
pub type MMHashSet<K> = HashSet<K, MMBuildHasher>;

// Take a bit-encoded k-mer (k <= 32) and decode it as a string of ACGT

pub fn decode_kmer64(kmer: Kmer64, k: u8) -> String {
    let mut seq = String::new();
    for i in 0..k {
        let c = (kmer >> (i * 2)) & 0b11;
        seq.push(match c {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        });
    }
    //reverse string
    seq.chars().rev().collect()
}

pub fn decode_kmer48(kmer: Kmer48, k: u8) -> String {
    let kmer = kmer.to_u64();
    let mut seq = String::new();
    for i in 0..k {
        let c = (kmer >> (i * 2)) & 0b11;
        seq.push(match c {
            0 => 'A',
            1 => 'C',
            2 => 'G',
            3 => 'T',
            _ => unreachable!(),
        });
    }
    //reverse string
    seq.chars().rev().collect()
}

#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct PreFragment {
    pub kmers_with_refpos: Vec<(Kmer64, u64)>,
    pub upper_base: usize,
    pub lower_base: usize,
    pub id: String,
}

#[derive(Debug, Default, Clone, Serialize, Deserialize, PartialEq, Eq)]
pub struct VarmerFragment {
    pub upper: usize,
    pub lower: usize,
    pub varmers: FxHashSet<usize>,
    pub upper_base: usize,
    pub lower_base: usize,
    pub id: String,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash)]
pub struct Varmer {
    pub kmer: Kmer64,
    pub count: u32,
    pub pos: u64,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
pub struct BasePileup {
    pub ref_pos: u64,
    pub ref_base: u8,
    pub base_freqs: [u32; 4],
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Eq, Hash, Default)]
pub struct TigdexOverlap {
    pub tig1: usize,
    pub tig2: usize,
    pub tig1_start: usize,
    pub tig1_end: usize,
    pub tig2_start: usize,
    pub tig2_end: usize,
    pub shared_tig: usize,
    pub variable_roots: usize,
    pub variable_tigs: usize,
    pub chain_reverse: bool,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct ChainInfo {
    pub chain: Vec<Anchor>,
    pub reverse: bool,
    pub score: i32,
    pub large_indel: bool,
}

pub type EdgeIndex = usize;
pub type NodeIndex = usize;

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct TigRead {
    pub tig_seq: Vec<u32>,
    pub id: String,
}

pub type Percentage = f64;
pub type Fraction = f32;

#[derive(Debug, Clone, PartialEq, Default, Serialize, Deserialize)]
pub struct TwinRead {
    /// Minimizer (syncmer) start positions, encoded with `encode_positions`. Only positions
    /// are stored: the k-mers are recomputed from `dna_seq` by `minimizer_kmers()`.
    pub(crate) minimizer_positions_enc: Vec<u8>,
    /// SNPmer start positions (after high-frequency filtering), encoded the same way.
    pub(crate) snpmer_positions_enc: Vec<u8>,
    pub blockmer_positions: Vec<u32>,
    pub blockmer_canonical: Vec<bool>,
    pub id: String,
    pub base_id: String,
    pub k: u8,
    pub l: u8,
    pub dna_seq: Seq<Dna>,
    pub qual_seq: Option<Seq<QualCompact3>>,
    pub est_id: Option<Percentage>,
    /// 0-based index of the input file this read came from (for --pooled-samples)
    pub file_idx: u32,
}

/// Encodes ascending positions as LEB128 varints of the gaps between them (the first
/// relative to 0). Seeds are a few bases apart, so most positions take one byte instead of
/// four, and none of the k-mer bytes need to be stored at all.
pub(crate) fn encode_positions(positions: &[u32]) -> Vec<u8> {
    let mut encoded = Vec::with_capacity(positions.len() + positions.len() / 8);
    let mut prev = 0;
    for &pos in positions {
        debug_assert!(pos >= prev, "positions must be ascending");
        let mut delta = pos - prev;
        while delta >= 0x80 {
            encoded.push((delta as u8) | 0x80);
            delta >>= 7;
        }
        encoded.push(delta as u8);
        prev = pos;
    }
    encoded.shrink_to_fit();
    encoded
}

/// Iterates over positions produced by `encode_positions`.
pub(crate) fn decode_positions(encoded: &[u8]) -> impl Iterator<Item = u32> + '_ {
    let mut bytes = encoded.iter();
    let mut pos = 0u32;
    std::iter::from_fn(move || {
        let mut delta = 0u32;
        let mut shift = 0;
        loop {
            let byte = *bytes.next()?;
            delta |= ((byte & 0x7F) as u32) << shift;
            if byte & 0x80 == 0 {
                break;
            }
            shift += 7;
        }
        pos += delta;
        Some(pos)
    })
}

/// Every varint ends in exactly one byte without the continuation bit.
fn count_positions(encoded: &[u8]) -> usize {
    encoded.iter().filter(|&&byte| byte & 0x80 == 0).count()
}

// `raw_kmer` reads the 2-bit bases straight out of bio-seq's `usize` storage words.
const _: () = assert!(usize::BITS == 64, "savont requires a 64-bit target");

/// The `k` bases of `seq` starting at `pos`, packed two bits per base with the first base
/// in the lowest bits (bio-seq's layout: A=0, C=1, G=2, T=3).
///
/// This replaces `Kmer::unsafe_from_seqslice`, which goes through bitvec's generic
/// `load_le` and was ~3x slower. That mattered once k-mers stopped being stored and every
/// lookup became a recompute.
#[inline]
fn raw_kmer(seq: &Seq<Dna>, pos: usize, k: usize) -> u64 {
    debug_assert!(k <= 32 && pos + k <= seq.len());
    let words = seq.into_raw();
    let bit = 2 * pos;
    let (word, offset) = (bit / 64, bit % 64);
    let mut raw = (words[word] as u64) >> offset;
    if offset + 2 * k > 64 {
        raw |= (words[word + 1] as u64) << (64 - offset);
    }
    let raw = raw & (u64::MAX >> (64 - 2 * k));
    // Assumes the bit vector starts at bit 0 of its first word, true for every `Seq` built
    // from bytes or cloned. Checked against bio-seq's own slice loading in debug builds.
    debug_assert_eq!(Ok(raw as usize), usize::try_from(&seq[pos..pos + k]));
    raw
}

#[derive(Clone, Copy, Debug, PartialEq, Eq, PartialOrd, Ord, Hash)]
#[repr(u8)]
pub enum QualCompact3 {
    Q33 = 0b0000,
    Q36 = 0b0001,
    Q39 = 0b0010,
    Q42 = 0b0011,
    Q45 = 0b0100,
    Q48 = 0b0101,
    Q51 = 0b0110,
    Q54 = 0b0111,
    Q57 = 0b1000,
    Q60 = 0b1001,
    Q63 = 0b1010,
    Q66 = 0b1011,
    Q69 = 0b1100,
    Q72 = 0b1101,
    Q75 = 0b1110,
    Q78 = 0b1111,
}

impl Codec for QualCompact3 {
    const BITS: u8 = 4;

    /// Take the two least significant bits of a `u8` and map them to the
    /// corresponding nucleotides.
    fn unsafe_from_bits(b: u8) -> Self {
        unsafe { std::mem::transmute(b & 0b1111) }
    }

    /// We can efficient verify that a byte is a valid `Dna` value if it's
    /// between 0 and 3.
    fn try_from_bits(b: u8) -> Option<Self> {
        // Round to nearest 3 and map to enum variant
        let rounded = match b {
            0..=34 => 0,   // Q33
            35..=37 => 1,  // Q36
            38..=40 => 2,  // Q39
            41..=43 => 3,  // Q42
            44..=46 => 4,  // Q45
            47..=49 => 5,  // Q48
            50..=52 => 6,  // Q51
            53..=55 => 7,  // Q54
            56..=58 => 8,  // Q57
            59..=61 => 9,  // Q60
            62..=64 => 10, // Q63
            65..=67 => 11, // Q66
            68..=70 => 12, // Q69
            71..=73 => 13, // Q72
            74..=76 => 14, // Q75
            _ => 15,       // Q78 or higher
        };

        // Use match instead of transmute for safety
        let m = match rounded {
            0 => Some(Self::Q33),
            1 => Some(Self::Q36),
            2 => Some(Self::Q39),
            3 => Some(Self::Q42),
            4 => Some(Self::Q45),
            5 => Some(Self::Q48),
            6 => Some(Self::Q51),
            7 => Some(Self::Q54),
            8 => Some(Self::Q57),
            9 => Some(Self::Q60),
            10 => Some(Self::Q63),
            11 => Some(Self::Q66),
            12 => Some(Self::Q69),
            13 => Some(Self::Q72),
            14 => Some(Self::Q75),
            15 => Some(Self::Q78),
            _ => None, // This case should never happen given our match above
        };

        return m;
    }

    /// The ASCII values of 'A', 'C', 'G', and 'T' can be translated into
    /// the numbers 0, 1, 2, and 3 using bitwise operations: `((b << 1) + b) >> 3`.
    fn unsafe_from_ascii(b: u8) -> Self {
        Self::unsafe_from_bits(b)
    }

    fn try_from_ascii(c: u8) -> Option<Self> {
        Self::try_from_bits(c)
    }

    //Not -33
    fn to_char(self) -> char {
        match self {
            QualCompact3::Q33 => '!',
            QualCompact3::Q36 => '"',
            QualCompact3::Q39 => '#',
            QualCompact3::Q42 => '$',
            QualCompact3::Q45 => '%',
            QualCompact3::Q48 => '&',
            QualCompact3::Q51 => '\'',
            QualCompact3::Q54 => '(',
            QualCompact3::Q57 => ')',
            QualCompact3::Q60 => '*',
            QualCompact3::Q63 => '+',
            QualCompact3::Q66 => ',',
            QualCompact3::Q69 => '-',
            QualCompact3::Q72 => '.',
            QualCompact3::Q75 => '/',
            QualCompact3::Q78 => '0',
        }
    }

    fn to_bits(self) -> u8 {
        self as u8
    }

    fn items() -> impl Iterator<Item = Self> {
        vec![
            QualCompact3::Q33,
            QualCompact3::Q36,
            QualCompact3::Q39,
            QualCompact3::Q42,
            QualCompact3::Q45,
            QualCompact3::Q48,
            QualCompact3::Q51,
            QualCompact3::Q54,
            QualCompact3::Q57,
            QualCompact3::Q60,
            QualCompact3::Q63,
            QualCompact3::Q66,
            QualCompact3::Q69,
            QualCompact3::Q72,
            QualCompact3::Q75,
            QualCompact3::Q78,
        ]
        .into_iter()
    }
}

impl ComplementMut for QualCompact3 {
    fn comp(&mut self) {}
}

impl Complement for QualCompact3 {}

#[inline]
fn reverse_bit_pairs(n: u64, k: usize) -> u64 {
    let even_mask: u64 = 0xAAAAAAAAAAAAAAAAu64;
    let odd_mask: u64 = 0x5555555555555555u64;

    let odd_bits_rev = (n & odd_mask).reverse_bits() >> (64 - 2 * k);
    let even_bits_rev = (n & even_mask).reverse_bits() >> (64 - 2 * k);

    return (odd_bits_rev >> 1) | (even_bits_rev << 1);
}

impl TwinRead {
    /// The k-mer at `pos` on the given strand (`canonical` = forward).
    #[inline]
    pub fn kmer_from_position_canonical(&self, pos: u32, k: usize, canonical: bool) -> Kmer48 {
        let raw = raw_kmer(&self.dna_seq, pos as usize, k);
        if canonical {
            Kmer48::from_u64(reverse_bit_pairs(raw, k))
        } else {
            Kmer48::from_u64(!raw & (u64::MAX >> (64 - 2 * k)))
        }
    }

    /// The k-mer at `pos` in seeding's canonical orientation: the strand that is smaller with
    /// the middle base masked out. Seeding skips k-mers whose two masked strands are equal
    /// (reverse-complement palindromic flanks), so every stored position has exactly one
    /// canonical k-mer and this returns the same k-mer seeding saw.
    #[inline]
    pub fn kmer_from_position(&self, pos: u32, k: usize) -> Kmer48 {
        let raw = raw_kmer(&self.dna_seq, pos as usize, k);
        // Our layout puts the first base in the highest bits, the reverse of bio-seq's...
        let forward = reverse_bit_pairs(raw, k);
        // ...so in our layout the reverse complement is simply bio-seq's bits complemented.
        let reverse = !raw & (u64::MAX >> (64 - 2 * k));
        let mid_mask = !(3 << (k - 1));
        if reverse & mid_mask < forward & mid_mask {
            Kmer48::from_u64(reverse)
        } else {
            Kmer48::from_u64(forward)
        }
    }

    pub fn minimizer_positions(&self) -> impl Iterator<Item = u32> + '_ {
        decode_positions(&self.minimizer_positions_enc)
    }

    pub fn snpmer_positions(&self) -> impl Iterator<Item = u32> + '_ {
        decode_positions(&self.snpmer_positions_enc)
    }

    pub fn minimizer_count(&self) -> usize {
        count_positions(&self.minimizer_positions_enc)
    }

    pub fn snpmer_count(&self) -> usize {
        count_positions(&self.snpmer_positions_enc)
    }

    pub(crate) fn set_minimizer_positions(&mut self, positions: &[u32]) {
        self.minimizer_positions_enc = encode_positions(positions);
    }

    pub(crate) fn set_snpmer_positions(&mut self, positions: &[u32]) {
        self.snpmer_positions_enc = encode_positions(positions);
    }

    /// Recomputes the minimizer k-mers from the sequence. Every call redoes the work, so a
    /// loop that visits the same read more than once should call this once and keep the result.
    pub fn minimizer_kmers(&self) -> Vec<Kmer48> {
        let k = self.k as usize;
        self.minimizer_positions()
            .map(|pos| self.kmer_from_position(pos, k))
            .collect()
    }

    /// Recomputes the SNPmer k-mers; the same caching advice as `minimizer_kmers` applies.
    pub fn snpmer_kmers(&self) -> Vec<Kmer48> {
        self.snpmers().map(|(_, kmer)| kmer).collect()
    }

    pub(crate) fn snpmers(&self) -> impl Iterator<Item = (u32, Kmer48)> + '_ {
        let k = self.k as usize;
        self.snpmer_positions()
            .map(move |pos| (pos, self.kmer_from_position(pos, k)))
    }

    pub fn snpmers_vec(&self) -> Vec<(u32, Kmer48)> {
        self.snpmers().collect()
    }

    /// Keeps only the SNPmers whose k-mer passes `keep`.
    pub fn retain_snpmers(&mut self, mut keep: impl FnMut(Kmer48) -> bool) {
        let kept: Vec<u32> = self
            .snpmers()
            .filter(|&(_, kmer)| keep(kmer))
            .map(|(pos, _)| pos)
            .collect();
        self.set_snpmer_positions(&kept);
    }

    /// Drops growth slack from the per-read buffers.
    ///
    /// bio-seq builds a `Seq` from a `Vec<u8>` without a size hint, so its buffer doubles
    /// past the data (~520 B for ~340 B of bases). Cloning a `Seq` allocates exactly. Twin
    /// reads live for the whole run, so one extra allocation per buffer at construction is
    /// worth it.
    pub fn shrink_to_fit(&mut self) {
        self.minimizer_positions_enc.shrink_to_fit();
        self.snpmer_positions_enc.shrink_to_fit();
        self.blockmer_positions.shrink_to_fit();
        self.blockmer_canonical.shrink_to_fit();
        self.dna_seq = self.dna_seq.clone();
        if let Some(qual) = &self.qual_seq {
            self.qual_seq = Some(qual.clone());
        }
    }

    pub fn blockmers_vec(&self) -> Vec<(u32, u64)> {
        self.blockmers().collect()
    }

    pub(crate) fn blockmers(&self) -> impl Iterator<Item = (u32, u64)> + '_ {
        self.blockmer_positions
            .iter()
            .zip(self.blockmer_canonical.iter())
            .map(|(&pos, &is_canonical)| {
                let kmer = self.kmer_from_position_canonical(
                    pos,
                    self.k as usize + self.l as usize,
                    is_canonical,
                );
                (pos, kmer.to_u64())
            })
    }
}

/// LSH bucket signatures of a read's minimizers, one per hash table.
pub(crate) fn lsh_signatures(minimizers: &[Kmer48]) -> [Option<u64>; LSH_NUM_TABLES] {
    use fxhash::FxHasher64;
    use std::hash::Hasher;

    if minimizers.len() < LSH_BUCKET_SIZE {
        return [None; LSH_NUM_TABLES];
    }
    std::array::from_fn(|table_idx| {
        // The old stable sort broke hash ties by the original minimizer order.
        let mut lowest = [(u64::MAX, usize::MAX); LSH_BUCKET_SIZE];
        for (index, kmer) in minimizers.iter().enumerate() {
            let mut hasher = FxHasher64::default();
            hasher.write_u64(table_idx as u64);
            hasher.write_u64(kmer.to_u64());
            let entry = (hasher.finish(), index);
            let mut slot = LSH_BUCKET_SIZE - 1;
            if entry < lowest[slot] {
                while slot > 0 && entry < lowest[slot - 1] {
                    lowest[slot] = lowest[slot - 1];
                    slot -= 1;
                }
                lowest[slot] = entry;
            }
        }
        let mut signature = 0;
        for (rank, (_, index)) in lowest.into_iter().enumerate() {
            signature ^= minimizers[index].to_u64().wrapping_mul(rank as u64 + 1);
        }
        Some(signature)
    })
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct KmerGlobalInfo {
    pub snpmer_info: Vec<SnpmerInfo>,
    pub solid_kmers: HashSet<Kmer48>,
    pub use_solid_kmers: bool,
    pub high_freq_thresh: f64,
    pub high_freq_kmers: HashSet<Kmer48>,
    pub read_files: Vec<PathBuf>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct TwinReadContainer {
    pub twin_reads: Vec<TwinRead>,
    pub outer_indices: Vec<usize>,
    //Not implemented TODO
    pub tig_reads: Vec<TwinRead>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Eq, Ord, PartialOrd, Hash)]
pub struct SnpmerInfo {
    pub split_kmer: u64,
    pub mid_bases: SmallVec<[u8; 2]>,
    pub counts: SmallVec<[u32; 2]>,
    pub k: u8,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Eq, Ord, PartialOrd, Hash)]
pub struct BlockmerInfo {
    /// The anchor k-mer (first k bases, without suffix)
    pub anchor_kmer: u64,
    /// The two most abundant blockmers (biallelic, with full sequence)
    pub blockmers: SmallVec<[Blockmer; 2]>,
    /// Counts for each blockmer
    pub counts: SmallVec<[u32; 2]>,
    /// Anchor k-mer length
    pub k: u8,
    /// Suffix length
    pub l: u8,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default)]
pub struct BlockmerGlobalInfo {
    pub blockmer_info: Vec<BlockmerInfo>,
    pub read_files: Vec<PathBuf>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub struct TwinOverlap {
    pub i1: usize,
    pub i2: usize,
    pub start1: usize,
    pub end1: usize,
    pub start2: usize,
    pub end2: usize,
    pub shared_minimizers: usize,
    pub shared_snpmers: usize,
    pub snpmers_in_both: (usize, usize),
    pub diff_snpmers: usize,
    pub chain_reverse: bool,
    pub chain_score: i32,
    pub large_indel: bool,
    pub minimizer_chain: Option<Vec<Anchor>>,
}

#[derive(Debug, Clone, PartialEq, Eq, Hash, Serialize, Deserialize, Default)]
pub struct SnpmerHit {
    pub pos1: u32,
    pub pos2: u32,
    pub bases: (u8, u8),
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Hash, Eq)]
pub struct CountsAndBases {
    pub counts: SmallVec<[[u32; 2]; 2]>,
    pub bases: SmallVec<[u8; 4]>,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Eq, Hash)]
pub struct AnchorBuilder {
    pub i: u32,
    pub j: u32,
    pub pos1: u32,
}

#[derive(Debug, Clone, Serialize, Deserialize, PartialEq, Default, Eq, Hash)]
pub struct Anchor {
    pub i: u32,
    pub j: u32,
    pub pos1: u32,
    pub pos2: u32,
}

// Enum for marking the state of a node during processing
#[derive(PartialEq, Eq, Clone, Debug, Hash, Copy)]
pub enum Direction {
    Incoming,
    Outgoing,
}
impl Direction {
    pub fn reverse(&self) -> Direction {
        match self {
            Direction::Incoming => Direction::Outgoing,
            Direction::Outgoing => Direction::Incoming,
        }
    }
}

#[inline]
pub fn bits_to_ascii(bit_rep: u8) -> u8 {
    match bit_rep {
        0 => b'A',
        1 => b'C',
        2 => b'G',
        3 => b'T',
        _ => unreachable!(),
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BareMappingOverlap {
    pub snpmer_identity: Fraction,
}

impl Eq for BareMappingOverlap {}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct TwoCycle {
    pub read_i: usize,
    pub read_j: usize,
    pub hang_penalty: i64,
    pub circular_length: usize,
    pub total_mini: usize,
}

#[derive(Debug, Clone, PartialEq, Default, Eq)]
pub struct BareInterval {
    pub start: u32,
    pub stop: u32,
}

impl Ord for BareInterval {
    #[inline]
    fn cmp(&self, other: &BareInterval) -> Ordering {
        match self.start.cmp(&other.start) {
            Ordering::Less => Ordering::Less,
            Ordering::Greater => Ordering::Greater,
            Ordering::Equal => self.stop.cmp(&other.stop),
        }
    }
}

impl PartialOrd for BareInterval {
    #[inline]
    fn partial_cmp(&self, other: &BareInterval) -> Option<Ordering> {
        Some(self.cmp(other))
    }
}

#[derive(Debug, Clone, PartialEq)]
pub struct Breakpoints {
    pub pos1: usize,
    pub pos2: usize,
    pub cov: usize,
    pub condition: i64,
}
#[derive(Debug, Clone, PartialEq, Default)]
pub struct GetSequenceInfoConfig {
    pub blunted: bool,
    pub dna_seq_info: bool,
}

#[derive(Debug, Clone, PartialEq)]
pub struct BubblePopResult {
    pub original_direction: Direction,
    pub end_direction: Direction,
    pub source_hash_id: NodeIndex,
    pub sink_hash_id: NodeIndex,
    pub remove_nodes: Vec<NodeIndex>,
    pub remove_edges: FxHashSet<EdgeIndex>,
}

impl BubblePopResult {
    pub fn new(
        original_direction: Direction,
        end_direction: Direction,
        source_hash_id: NodeIndex,
        sink_hash_id: NodeIndex,
        remove_nodes: Vec<NodeIndex>,
        remove_edges: FxHashSet<EdgeIndex>,
    ) -> Self {
        BubblePopResult {
            original_direction,
            end_direction,
            source_hash_id,
            sink_hash_id,
            remove_nodes,
            remove_edges,
        }
    }
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct BeamSearchSoln {
    pub path: Vec<EdgeIndex>,
    pub coverages: Vec<(MultiCov, usize)>,
    pub score: f64,
    pub path_nodes: Vec<NodeIndex>,
    pub depth: usize,
    pub current_length: usize,
}

pub struct BeamStartState {
    pub initial_unitig_length: usize,
    pub initial_unitig_size: usize,
}

#[derive(Debug, Clone, PartialEq, Default)]
pub struct OverlapAdjMap {
    pub adj_map: FxHashMap<NodeIndex, Vec<NodeIndex>>,
}

pub fn dna_seq_to_u8(slice: &Seq<Dna>) -> Vec<u8> {
    slice
        .iter()
        .map(|x| x.to_char().to_ascii_uppercase() as u8)
        .collect()
}

pub fn dna_slice_to_u8(slice: &SeqSlice<Dna>) -> Vec<u8> {
    slice
        .iter()
        .map(|x| x.to_char().to_ascii_uppercase() as u8)
        .collect()
}

pub fn quality_slice_to_u8(slice: &SeqSlice<QualCompact3>) -> Vec<u8> {
    slice.iter().map(|x| x as u8 * 3).collect()
}

pub fn quality_seq_to_u8(slice: &Seq<QualCompact3>) -> Vec<u8> {
    slice.iter().map(|x| x as u8 * 3).collect()
}

pub fn revcomp_u8(seq: &Vec<u8>) -> Vec<u8> {
    seq.iter()
        .rev()
        .map(|x| match x {
            b'A' => b'T',
            b'T' => b'A',
            b'C' => b'G',
            b'G' => b'C',
            _ => b'N',
        })
        .collect()
}

#[derive(Debug, Clone, PartialEq)]
pub struct CompareTwinReadOptions {
    pub compare_snpmers: bool,
    pub retain_chain: bool,
    pub force_query_nonoverlap: bool,
    pub force_ref_nonoverlap: bool,
    pub supplementary_threshold_score: Option<f64>,
    pub supplementary_threshold_ratio: Option<f64>,
    // When not forcing 1-to-1 alignments, allow query overlaps only if secondary threshold is below a certain amount
    pub secondary_threshold: Option<f64>,
    //Preload
    pub read1_mininimizers: Option<Vec<(u32, Kmer48)>>,
    pub read1_snpmers: Option<Vec<(u32, Kmer48)>>,
    pub max_gap: usize,
    pub double_gap: usize,
}

impl Default for CompareTwinReadOptions {
    fn default() -> Self {
        CompareTwinReadOptions {
            compare_snpmers: true,
            retain_chain: false,
            force_query_nonoverlap: false,
            force_ref_nonoverlap: true,
            supplementary_threshold_score: Some(500.0),
            supplementary_threshold_ratio: Some(0.25),
            secondary_threshold: Some(0.50),
            read1_mininimizers: None,
            read1_snpmers: None,
            max_gap: MAX_GAP_CHAINING,
            double_gap: 10_000,
        }
    }
}

pub struct HeavyCutOptions<'a> {
    pub samples: usize,
    pub temperature: f64,
    pub steps: usize,
    pub max_forward: usize,
    pub max_reads_forward: usize,
    pub safe_length_back: usize,
    pub ol_thresh: f64,
    pub tip_threshold: usize,
    pub strain_repeat_map: Option<&'a FxHashMap<NodeIndex, FxHashSet<NodeIndex>>>,
    pub special_small: bool,
    pub max_length_search: usize,
    pub require_safety: bool,
    pub only_tips: bool,
    pub cut_tips: bool,
    pub debug: bool,
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn lsh_signatures_match_stable_sort_with_repeated_minimizers() {
        use fxhash::FxHasher64;
        use std::hash::Hasher;

        let raw = [9, 2, 4, 2, 1, 9, 38, 1, 14, 15, 64, 38, 0xFFFF_FFFF_FFFF];
        for end in 0..=raw.len() {
            let minimizers: Vec<Kmer48> =
                raw[..end].iter().copied().map(Kmer48::from_u64).collect();
            let expected = std::array::from_fn(|table_idx| {
                if end < LSH_BUCKET_SIZE {
                    return None;
                }
                let mut hashed: Vec<_> = raw[..end]
                    .iter()
                    .copied()
                    .map(|kmer| {
                        let mut hasher = FxHasher64::default();
                        hasher.write_u64(table_idx as u64);
                        hasher.write_u64(kmer);
                        (hasher.finish(), kmer)
                    })
                    .collect();
                hashed.sort_by_key(|&(hash, _)| hash);
                let signature = hashed[..LSH_BUCKET_SIZE]
                    .iter()
                    .enumerate()
                    .fold(0, |value, (i, &(_, kmer))| {
                        value ^ kmer.wrapping_mul(i as u64 + 1)
                    });
                Some(signature)
            });
            assert_eq!(lsh_signatures(&minimizers), expected);
        }
    }

    #[test]
    fn positions_round_trip() {
        let positions = [0, 1, 5, 127, 128, 300, 16_683, 16_684, 2_000_000, u32::MAX];
        let encoded = encode_positions(&positions);
        assert_eq!(decode_positions(&encoded).collect::<Vec<_>>(), positions);
        assert_eq!(count_positions(&encoded), positions.len());
        assert!(encode_positions(&[]).is_empty());
        assert_eq!(decode_positions(&[]).count(), 0);
    }

    /// Recomputed k-mers must equal what seeding computes from the raw bytes: forward strand
    /// when its masked form is strictly smaller, otherwise the reverse complement.
    #[test]
    fn kmer_from_position_matches_seeding() {
        let mut state = 0x9E37_79B9_7F4A_7C15u64;
        let bases: Vec<u8> = (0..1000)
            .map(|_| {
                state ^= state << 13;
                state ^= state >> 7;
                state ^= state << 17;
                b"ACGT"[(state >> 32) as usize % 4]
            })
            .collect();
        let read = TwinRead {
            dna_seq: bases.clone().try_into().unwrap(),
            ..TwinRead::default()
        };
        for k in [13, 15, 17, 19, 21, 23] {
            let mask = u64::MAX >> (64 - 2 * k);
            let mid_mask = !(3u64 << (k - 1));
            for pos in 0..=bases.len() - k {
                let forward = bases[pos..pos + k]
                    .iter()
                    .fold(0u64, |acc, &b| (acc << 2) | BYTE_TO_SEQ[b as usize] as u64);
                let reverse = bases[pos..pos + k]
                    .iter()
                    .rev()
                    .fold(0u64, |acc, &b| (acc << 2) | (3 - BYTE_TO_SEQ[b as usize] as u64))
                    & mask;
                let pos32 = pos as u32;
                assert_eq!(read.kmer_from_position_canonical(pos32, k, true).to_u64(), forward);
                assert_eq!(read.kmer_from_position_canonical(pos32, k, false).to_u64(), reverse);
                if forward & mid_mask == reverse & mid_mask {
                    continue; // palindromic flanks: seeding never stores these
                }
                let seeding = if forward & mid_mask < reverse & mid_mask { forward } else { reverse };
                assert_eq!(read.kmer_from_position(pos32, k).to_u64(), seeding);
            }
        }
    }

    #[test]
    fn convert_from_u64() {
        let kmer = Kmer48::from_u64(1);
        assert_eq!(kmer.data, [1, 0, 0, 0, 0, 0]);
    }

    #[test]
    fn bioseq_vs_ours() {
        let kmer_bioseq = Kmer::<Dna, 5, u64>::unsafe_from_seqslice(dna!("ACGTG"));
        let kmer_bioseq_u64 = kmer_bioseq.bs;
        println!("{:08b}", kmer_bioseq_u64);
        let kmer_ours = 0b00_01_10_11_10;
        let reversed_bioseq = reverse_bit_pairs(kmer_bioseq_u64, 5);
        println!("{:08b}", reversed_bioseq);

        assert_eq!(reverse_bit_pairs(kmer_bioseq_u64, 5), kmer_ours);
    }

    #[test]
    fn reverse_comp_kmer() {
        // ACGTG
        let kmer_ours = 0b00_01_10_11_10;
        //reverse comp is CACGT
        let goal = 0b01_00_01_10_11;
        let kmer_ours_rev_comp = reverse_bit_pairs(kmer_ours ^ (u64::MAX), 5);

        assert_eq!(kmer_ours_rev_comp, goal);
    }
}
