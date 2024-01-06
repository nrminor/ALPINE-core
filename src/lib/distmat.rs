use clap::ValueEnum;
use displaydoc::Display;
use std::io;
use textdistance::{
    nstr::{lcsseq, lcsstr},
    str::{damerau_levenshtein, jaro_winkler, levenshtein, ratcliff_obershelp, smith_waterman},
    str::{hamming, jaccard},
};

#[derive(ValueEnum, Debug, Clone, PartialEq, Display)]
pub enum DistanceMethods {
    /// Hamming edit distance
    Hamming,

    /// Levenshtein edit distance
    Levenshtein,

    /// Damerau-Levenshtein edit distance
    DamerauLevenshtein,

    /// Jaro-Winkler edit distance
    JaroWinkler,

    /// Smith-Waterman edit distance
    SmithWaterman,

    /// Ratcliff-Obershelp/Gestalt pattern matching sequence-based distance
    RatcliffObershelp,

    /// Longest Common SubSequence distance
    LCSSeq,

    /// Longest Common SubString
    LCSStr,

    /// Jaccard token/kmer-based distance
    Jaccard,
}

#[derive(ValueEnum, Debug, Clone, PartialEq)]
pub enum Stringency {
    Lenient,
    Intermediate,
    Strict,
    Extreme,
}

trait DistanceCalculator {
    fn calculate_distance(&self, s1: &str, s2: &str) -> f64;
}

impl DistanceCalculator for DistanceMethods {
    fn calculate_distance(&self, s1: &str, s2: &str) -> f64 {
        match self {
            DistanceMethods::Hamming => hamming(s1, s2) as f64,
            DistanceMethods::Levenshtein => levenshtein(s1, s2) as f64,
            DistanceMethods::DamerauLevenshtein => damerau_levenshtein(s1, s2) as f64,
            DistanceMethods::JaroWinkler => jaro_winkler(s1, s2),
            DistanceMethods::SmithWaterman => smith_waterman(s1, s2) as f64,
            DistanceMethods::RatcliffObershelp => ratcliff_obershelp(s1, s2),
            DistanceMethods::LCSSeq => lcsseq(s1, s2),
            DistanceMethods::LCSStr => lcsstr(s1, s2),
            DistanceMethods::Jaccard => jaccard(s1, s2),
        }
    }
}

fn _weight_by_cluster_size() {}

pub fn compute_distance_matrix() -> io::Result<()> {
    println!("API for computing a pairwise distance matrix using one of a few kinds of edit distances coming soon!");
    Ok(())
}
