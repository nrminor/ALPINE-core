use anyhow::Result;
use clap::ValueEnum;
use displaydoc::Display;
use polars::{lazy::dsl::col, prelude::*};
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

fn _weight_by_cluster_size(
    seq_name: &str,
    stringency: &str,
    dist_df: &LazyFrame,
    cluster_table: &LazyFrame,
) -> Result<LazyFrame> {
    // separate out the columns of information we need
    let cluster_query = cluster_table.clone().collect()?;
    let col_names = &cluster_query.get_column_names();
    let member_types = col_names.first().unwrap();
    let seq_names = col_names.get(8).unwrap();
    let cluster_sizes = col_names.get(2).unwrap();

    // Filter down the df so that only rows representing centroids are present
    let centroid_lf = cluster_table
        .clone()
        .filter(col(member_types).eq(lit("C")))
        .sort(seq_name, Default::default());

    // extract a frame of cluster sizes that are in the same order as the sequences in the distance matrix
    let all_sizes = centroid_lf
        .clone()
        .select([col(cluster_sizes), col(seq_names)]);

    // Filter down to hits only and use to get a total number of sequences for the current month
    let hits_df = cluster_table
        .clone()
        .filter(col(member_types).eq(lit("H")))
        .sort(seq_name, Default::default())
        .collect()?;
    let month_total = if hits_df.shape().0 == 0 {
        1.0
    } else {
        hits_df.shape().0 as f64
    };

    // find cluster size for the accession in question
    let weighting_size = centroid_lf
        .clone()
        .filter(col(seq_names).eq(lit(seq_name)))
        .select([col(cluster_sizes)])
        .collect()?
        .column(cluster_sizes)?
        .iter()
        .next()
        .unwrap()
        .try_extract::<f64>()?;
    let weighting_freq = weighting_size / month_total;

    // compute weights lazyframe column
    let weights_lf = if stringency == "strict" || stringency == "extreme" {
        all_sizes.with_column(
            (col(cluster_sizes) * lit(weighting_freq.ln() * -1.0))
                / lit(month_total).alias("Weights"),
        )
    } else {
        all_sizes.with_column((lit(1) - col(cluster_sizes)) / lit(month_total).alias("Weights"))
    };

    // double check that we now have as many weights as we need
    assert!(
        weights_lf.clone().collect()?.shape().0 == dist_df.clone().collect()?.shape().0,
        "ALPINE was not able to find the same number of clusters as are represented in the centroid distance matrix"
    );

    Ok(weights_lf)
}

pub fn compute_distance_matrix(
    _fasta: &str,
    _cluster_table: &str,
    _yearmonth: &str,
    _stringency: &Stringency,
    _distance_method: &DistanceMethods,
) -> io::Result<()> {
    println!("API for computing a pairwise distance matrix using one of a few kinds of edit distances coming soon!");
    // DistMatrix::from_pw_distances_with(&fasta_seq_vec, |seq1, seq2| distance_method.calculate_distance(seq1, seq2));
    Ok(())
}
