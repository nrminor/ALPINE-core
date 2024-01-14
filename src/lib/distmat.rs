use anyhow::Result;
use clap::ValueEnum;
use displaydoc::Display;
use distmat::SquareMatrix;
use noodles::{bgzf, fasta};
use polars::{lazy::dsl::col, prelude::*};
use std::fs::File;
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

fn weight_by_cluster_size(
    seq_name: &str,
    stringency: &Stringency,
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
    let weights_lf = match *stringency {
        Stringency::Strict | Stringency::Extreme => all_sizes.with_column(
            (col(cluster_sizes) * lit(weighting_freq.ln() * -1.0))
                / lit(month_total).alias("Weights"),
        ),
        _ => {
            all_sizes.with_column((lit(1) - col(cluster_sizes)) / lit(month_total).alias("Weights"))
        }
    };

    // double check that we now have as many weights as we need
    assert!(
        weights_lf.clone().collect()?.shape().0 == dist_df.clone().collect()?.shape().0,
        "ALPINE was not able to find the same number of clusters as are represented in the centroid distance matrix"
    );

    Ok(weights_lf)
}

pub fn compute_distance_matrix(
    fasta: &str,
    cluster_table: &str,
    yearmonth: &str,
    stringency: &Stringency,
    distance_method: &DistanceMethods,
) -> Result<()> {
    println!("API for computing a pairwise distance matrix using one of a few kinds of edit distances coming soon!");

    // buffer the incoming fasta file
    let mut fa_reader = File::open(fasta)
        .map(bgzf::Reader::new)
        .map(fasta::Reader::new)?;

    // collect the FASTA IDs and sequences
    let (ids, sequences): (Vec<String>, Vec<String>) = fa_reader
        .records()
        .map(|result| {
            let record = result.expect("Error during fasta record parsing.");
            let sequence_string = String::from_utf8(record.sequence().get(..).unwrap().to_vec())
                .expect("Found invalid UTF-8 in sequence");
            (record.name().to_owned(), sequence_string)
        })
        .unzip();

    // call a distance matrix with the chosen distance metric (defaulting to Levenshtein)
    let mut pw_distmat = SquareMatrix::from_pw_distances_with(&sequences, |seq1, seq2| {
        distance_method.calculate_distance(seq1, seq2)
    });
    pw_distmat.set_labels(ids.clone());

    // pull distmat information out of the SquareMatrix struct and convert to dataframe
    let mut dist_col_vev: Vec<Series> = Vec::with_capacity(pw_distmat.size());
    for (i, (column, label)) in pw_distmat
        .iter_cols()
        .zip(pw_distmat.iter_labels())
        .enumerate()
    {
        let series = Series::new(label, column.collect::<Vec<f64>>());
        dist_col_vev[i] = series
    }
    let mut dist_lf = DataFrame::new(dist_col_vev)?.lazy();

    // read the cluster table into a lazyframe to query for cluster-size-based distance weights
    let cluster_df = CsvReader::from_path(cluster_table)?
        .has_header(false)
        .finish()?
        .lazy();

    // multiply weights onto each column based on the sequence it represents
    for id in ids {
        let weights_lf = weight_by_cluster_size(&id, stringency, &dist_lf, &cluster_df)?;
        dist_lf = dist_lf.join(
            weights_lf,
            [col(&id)],
            [col(&id)],
            JoinArgs::new(JoinType::Cross),
        );
    }
    let mut final_df = dist_lf.collect()?;

    // write out the weighted distance matrix
    let out_name = format!("{}-distmat.csv", yearmonth);
    let out_handle = File::create(out_name).unwrap();
    CsvWriter::new(out_handle)
        .has_header(true)
        .finish(&mut final_df)
        .unwrap();

    Ok(())
}
