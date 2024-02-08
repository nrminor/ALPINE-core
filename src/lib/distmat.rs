use anyhow::anyhow;
use anyhow::Result;
use bio::alignment::distance::simd::{hamming, levenshtein};
use block_aligner::{scan_block::*, scores::*};
use clap::ValueEnum;
use derive_new::new;
use distmat::SquareMatrix;
use noodles::{bgzf, fasta};
use polars::{lazy::dsl::col, prelude::*};
use rayon::prelude::*;
use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::io::ErrorKind;
use std::ops::Mul;
use std::sync::Arc;
use textdistance::{
    str::{damerau_levenshtein, jaro_winkler, ratcliff_obershelp, smith_waterman},
    str::{entropy_ncd, jaccard},
};

#[derive(ValueEnum, Debug, Clone, PartialEq)]
pub enum DistanceMethods {
    /// Default distance based on block-alignment
    BlockAlignment,

    /// Hamming edit distance (NOTE: Only works with aligned sequences)
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

    /// Jaccard token/kmer-based distance
    Jaccard,

    /// Entropy normalized compression distance
    Entropy,
}

impl fmt::Display for DistanceMethods {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                DistanceMethods::BlockAlignment => "block-alignment",
                DistanceMethods::Hamming => "hamming",
                DistanceMethods::Levenshtein => "levenshtein",
                DistanceMethods::DamerauLevenshtein => "damerau-levenshtein",
                DistanceMethods::JaroWinkler => "jaro-winkler",
                DistanceMethods::SmithWaterman => "smith-waterman",
                DistanceMethods::RatcliffObershelp => "ratcliff-obershelp",
                DistanceMethods::Jaccard => "jaccard",
                DistanceMethods::Entropy => "entropy",
            }
        )
    }
}

#[derive(ValueEnum, Debug, Clone, PartialEq)]
pub enum Stringency {
    Lenient,
    Intermediate,
    Strict,
    Extreme,
}

impl fmt::Display for Stringency {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                Stringency::Lenient => "lenient",
                Stringency::Intermediate => "intermediate",
                Stringency::Strict => "strict",
                Stringency::Extreme => "extreme",
            }
        )
    }
}

fn generate_score_matrix() -> NucMatrix {
    // generate mutable scoring matrix to use for constraining
    // the characters that will be counted toward mismatches.
    // Only A's, T's, G's, C's, U's, and -'s will count.
    let mut scoring_matrix = NucMatrix::new_simple(0, -1);

    /*
    KEY FOR AMBIGUOUS BASES
    =======================
    N: Any nucleotide (A or C or G or T or U)
    R: Purine (A or G)
    Y: Pyrimidine (T or C)
    K: Keto (G or T)
    M: Amino (A or C)
    S: Strong interaction (3 H bonds) (G or C)
    W: Weak interaction (2 H bonds) (A or T)
    B: Not A (C or G or T)
    D: Not C (A or G or T)
    H: Not G (A or C or T)
    V: Not T or U (A or C or G)
    */

    // account for masked "any" bases, N
    scoring_matrix.set(b'N', b'A', 0);
    scoring_matrix.set(b'N', b'T', 0);
    scoring_matrix.set(b'N', b'G', 0);
    scoring_matrix.set(b'N', b'C', 0);

    // account for any purine bases, R
    scoring_matrix.set(b'R', b'A', 0);
    scoring_matrix.set(b'R', b'T', 0);
    scoring_matrix.set(b'R', b'G', 0);
    scoring_matrix.set(b'R', b'C', 0);

    // account for any Pyrimidine bases, Y
    scoring_matrix.set(b'Y', b'A', 0);
    scoring_matrix.set(b'Y', b'T', 0);
    scoring_matrix.set(b'Y', b'G', 0);
    scoring_matrix.set(b'Y', b'C', 0);

    // account for any Keto bases, K
    scoring_matrix.set(b'K', b'A', 0);
    scoring_matrix.set(b'K', b'T', 0);
    scoring_matrix.set(b'K', b'G', 0);
    scoring_matrix.set(b'K', b'C', 0);

    // account for any amino bases, M
    scoring_matrix.set(b'M', b'A', 0);
    scoring_matrix.set(b'M', b'T', 0);
    scoring_matrix.set(b'M', b'G', 0);
    scoring_matrix.set(b'M', b'C', 0);

    // account for any strong interaction bases, S
    scoring_matrix.set(b'S', b'A', 0);
    scoring_matrix.set(b'S', b'T', 0);
    scoring_matrix.set(b'S', b'G', 0);
    scoring_matrix.set(b'S', b'C', 0);

    // account for and weak interaction bases, W
    scoring_matrix.set(b'W', b'A', 0);
    scoring_matrix.set(b'W', b'T', 0);
    scoring_matrix.set(b'W', b'G', 0);
    scoring_matrix.set(b'W', b'C', 0);

    // account for any "not A" bases, B
    scoring_matrix.set(b'B', b'A', 0);
    scoring_matrix.set(b'B', b'T', 0);
    scoring_matrix.set(b'B', b'G', 0);
    scoring_matrix.set(b'B', b'C', 0);

    // account for and "not C" bases, D
    scoring_matrix.set(b'D', b'A', 0);
    scoring_matrix.set(b'D', b'T', 0);
    scoring_matrix.set(b'D', b'G', 0);
    scoring_matrix.set(b'D', b'C', 0);

    // account for any "not G" bases, H
    scoring_matrix.set(b'H', b'A', 0);
    scoring_matrix.set(b'H', b'T', 0);
    scoring_matrix.set(b'H', b'G', 0);
    scoring_matrix.set(b'H', b'C', 0);

    // account for any "not T or U" bases, V
    scoring_matrix.set(b'V', b'A', 0);
    scoring_matrix.set(b'V', b'T', 0);
    scoring_matrix.set(b'V', b'G', 0);
    scoring_matrix.set(b'V', b'C', 0);

    scoring_matrix
}

fn block_distance(alpha: &[u8], beta: &[u8]) -> f64 {
    // set up parameters for alignment
    // let (min_block_size, max_block_size) = if beta.len() <= alpha.len() {
    //     (percent_len(beta.len(), 0.01), percent_len(beta.len(), 0.1))
    // } else {
    //     (
    //         percent_len(alpha.len(), 0.01),
    //         percent_len(alpha.len(), 0.1),
    //     )
    // };
    let (min_block_size, max_block_size) = (32, 256);

    // pack nucleotide bytes into SIMD-able blocks
    let reference = PaddedBytes::from_bytes::<NucMatrix>(alpha, max_block_size);
    let query = PaddedBytes::from_bytes::<NucMatrix>(beta, max_block_size);

    // set up parameters for how matches, mismatches, and gaps are scored
    let scoring_matrix = generate_score_matrix();
    let gap_penalty = Gaps {
        open: -2,
        extend: -1,
    };

    // instantiate an aligner
    let mut aligner = Block::<true, false>::new(query.len(), reference.len(), max_block_size);
    aligner.align(
        &query,
        &reference,
        &scoring_matrix,
        gap_penalty,
        min_block_size..=max_block_size,
        0,
    );

    // return the final score
    -aligner.res().score as f64
}

trait DistanceCalculator {
    fn calculate_distance(&self, s1: &str, s2: &str) -> f64;
}

impl DistanceCalculator for DistanceMethods {
    fn calculate_distance(&self, s1: &str, s2: &str) -> f64 {
        match self {
            DistanceMethods::BlockAlignment => block_distance(s1.as_bytes(), s2.as_bytes()),
            DistanceMethods::Hamming => hamming(s1.as_bytes(), s2.as_bytes()) as f64,
            DistanceMethods::Levenshtein => levenshtein(s1.as_bytes(), s2.as_bytes()) as f64,
            DistanceMethods::DamerauLevenshtein => damerau_levenshtein(s1, s2) as f64,
            DistanceMethods::JaroWinkler => jaro_winkler(s1, s2),
            DistanceMethods::SmithWaterman => smith_waterman(s1, s2) as f64,
            DistanceMethods::RatcliffObershelp => ratcliff_obershelp(s1, s2),
            DistanceMethods::Jaccard => 1.0 - jaccard(s1, s2),
            DistanceMethods::Entropy => entropy_ncd(s1, s2),
        }
    }
}

fn collect_fa_data(fasta: &str) -> Result<(Vec<String>, Vec<Arc<str>>)> {
    // pull out record IDs and sequences into their own string vectors, while
    // handling potential bgzip compression
    let parsed_fasta: std::io::Result<Vec<(String, Arc<str>)>> = if fasta.ends_with(".gz") {
        File::open(fasta)
            .map(bgzf::Reader::new)
            .map(fasta::Reader::new)?
            .records()
            .map(|result| {
                result.and_then(|record| {
                    let id = record.name().to_owned();
                    unpack_sequence(&record).map(|sequence_string| (id, Arc::from(sequence_string)))
                })
            })
            .collect()
    } else {
        File::open(fasta)
            .map(BufReader::new)
            .map(fasta::Reader::new)?
            .records()
            .map(|result| {
                result.and_then(|record| {
                    let id = record.name().to_owned();
                    unpack_sequence(&record).map(|sequence_string| (id, Arc::from(sequence_string)))
                })
            })
            .collect()
    };

    let (ids, sequences) = match parsed_fasta {
        Ok(pairs) => pairs.into_par_iter().unzip(),
        Err(e) => return Err(e.into()),
    };

    Ok((ids, sequences))
}

/// Cluster columns contains the column names where the information ALPINE needs is stored
#[derive(new, Debug, Clone)]
struct ClusterColumns {
    type_col: Arc<str>,
    index_col: Arc<str>,
    id_col: Arc<str>,
    size_col: Arc<str>,
}

fn get_cluster_cols(cluster_table: &LazyFrame) -> Result<ClusterColumns> {
    // separate out the columns of information we need
    let cluster_query = cluster_table.clone().collect()?;
    let col_names = cluster_query.get_column_names();

    let type_col: Arc<str> = match col_names.first() {
        Some(col_name) => Arc::from(col_name.to_string()),
        None => {
            eprintln!(
                "Please double check that the column of VSEARCH cluster types is the first column."
            );
            return Err(anyhow!(
                "Member types could not be parsed from provided cluster table,"
            ));
        }
    };

    let index_col: Arc<str> = match col_names.get(1) {
        Some(col_name) => Arc::from(col_name.to_string()),
        None => {
            eprintln!(
                "Please double check that the column of VSEARCH cluster index is the second column."
            );
            return Err(anyhow!(
                "Column indices could not be parsed from provided cluster table,"
            ));
        }
    };

    let name_col: Arc<str> = match col_names.get(8) {
        Some(col_name) => Arc::from(col_name.to_string()),
        None => {
            eprintln!("Please double check that the column of sequence names is the ninth column.");
            return Err(anyhow!(
                "Sequence names could not be parsed from provided cluster table,"
            ));
        }
    };

    let size_col: Arc<str> = match col_names.get(2) {
        Some(col_name) => Arc::from(col_name.to_string()),
        None => {
            eprintln!(
                "Please double check that the column of VSEARCH cluster sizes is the third column."
            );
            return Err(anyhow!(
                "Cluster sizes could not be parsed from provided cluster table,"
            ));
        }
    };

    Ok(ClusterColumns::new(type_col, index_col, name_col, size_col))
}

fn get_size_per_member(
    cluster_table: &LazyFrame,
    centroids_only: &LazyFrame,
    clust_cols: &ClusterColumns,
) -> Result<(f64, DataFrame)> {
    // pull out the sizes for each centroid by index
    let centroid_sizes = centroids_only
        .clone()
        .select(&[col(&clust_cols.index_col), col(&clust_cols.size_col)])
        .collect()?;

    // Filter down to hits only and use to get a total number of sequences for the current month
    let member_count = cluster_table
        .clone()
        .filter(
            col(&clust_cols.type_col)
                .eq(lit("H"))
                .or(col(&clust_cols.type_col).eq(lit("S"))),
        )
        .select(&[col(&clust_cols.index_col)])
        .collect()?
        .shape()
        .0;

    // count rows of members to get the month total
    let month_total: f64 = if member_count == 0 {
        1.0
    } else {
        member_count as f64
    };

    Ok((month_total, centroid_sizes))
}

fn get_cluster_index(
    cluster_table: &LazyFrame,
    clust_cols: &ClusterColumns,
    seq_name: &str,
) -> Result<i64> {
    let filtered = cluster_table
        .clone()
        .filter(col(&clust_cols.id_col).eq(lit(seq_name)))
        .select(&[col(&clust_cols.index_col)])
        .collect()?;

    let index = filtered
        .column(&clust_cols.index_col)?
        .get(0)?
        .try_extract::<i64>()?;

    Ok(index)
}

fn compute_weighting_freq(
    centroid_lf: LazyFrame,
    clust_index: i64,
    month_total: f64,
    clust_cols: &ClusterColumns,
) -> Result<f64> {
    let collected_df = centroid_lf
        .clone()
        .filter(col(&clust_cols.index_col).eq(clust_index))
        .select([col(&clust_cols.size_col)])
        .collect()?;

    let attempt = match collected_df
        .column(&clust_cols.size_col)?
        .iter()
        .next() {
            Some(value) => value,
            None => return Err(anyhow!("Could not parse centroid data to compute a weight. Please double check the input cluster table."))
        };

    let cluster_freq = attempt.try_extract::<f64>()? / month_total;

    Ok(cluster_freq)
}

fn weight_by_cluster_size(
    seq_name: &str,
    stringency: &Stringency,
    cluster_table: &LazyFrame,
) -> Result<(String, Series)> {
    let clust_cols = get_cluster_cols(cluster_table)?;

    // Filter down the df so that only rows representing centroids are present
    let centroids_only = cluster_table
        .clone()
        .filter(col(&clust_cols.type_col).eq(lit("C")));

    // Filter down to hits only and use to get a total number of sequences for the current month
    let (month_total, all_size_df) =
        get_size_per_member(cluster_table, &centroids_only, &clust_cols)?;

    // determine the cluster index for the current cluster member
    let index: i64 = get_cluster_index(cluster_table, &clust_cols, seq_name)?;

    // find the frequency of the cluster for the member accession in question
    let weighting_freq = compute_weighting_freq(centroids_only, index, month_total, &clust_cols)?;

    // compute Polars series of weights along with the series name
    let weights_header = format!("{}_weights", seq_name);
    let weights_lf = match *stringency {
        Stringency::Strict | Stringency::Extreme => all_size_df
            .lazy()
            .with_column(lit(-1.0).alias("negative"))
            .with_column(lit(weighting_freq.ln()).alias("log_freq"))
            .with_column(lit(month_total).alias("total"))
            .with_column(
                col(&clust_cols.size_col) * ((col("negative") * col("log_freq")) / col("total")),
            )
            .rename([&clust_cols.size_col], [&weights_header]),
        _ => all_size_df
            .lazy()
            .with_column(lit(1.0).alias("tmp_int"))
            .with_column(lit(-1.0).alias("negative"))
            .with_column(lit(weighting_freq).alias("freq"))
            .with_column(lit(month_total).alias("total"))
            .with_column(
                col(&clust_cols.size_col)
                    * ((col("tmp_int") + col("negative") * col("freq")) / col("total")),
            )
            .rename([&clust_cols.size_col], [&weights_header]),
    };

    let weights = weights_lf
        .select(&[col(&weights_header)])
        .collect()?
        .column(&weights_header)?
        .to_owned();

    Ok((weights_header, weights))
}

fn unpack_sequence(record: &fasta::Record) -> std::io::Result<String> {
    let seq_attempt =
        match record.sequence().get(..) {
            Some(seq) => seq.to_vec(),
            None => return Err(std::io::Error::new(
                ErrorKind::InvalidData,
                "No sequence was found for the provided record. Double check FASTA completeness.",
            )),
        };

    let seq_as_string = String::from_utf8(seq_attempt).unwrap();

    Ok(seq_as_string)
}

fn process_cluster_info(
    cluster_table: Option<&str>,
    dist_col_vec: Vec<Series>,
    ids: &Vec<String>,
    stringency: &Stringency,
) -> Result<DataFrame> {
    let mut dist_df = DataFrame::new(dist_col_vec)?;
    dist_df = match cluster_table {
        Some(table) => {
            // read the cluster table into a lazyframe to query for cluster-size-based distance weights
            let cluster_df = CsvReader::from_path(table)?
                .has_header(false)
                .with_delimiter(b'\t')
                .finish()?
                .lazy();

            // multiply weights onto each column based on the sequence it represents
            for id in ids {
                let (weights_header, weights) =
                    weight_by_cluster_size(id, stringency, &cluster_df)?;
                dist_df = dist_df
                    .hstack(&[weights])?
                    .lazy()
                    .with_columns(&[col(id).mul(col(&weights_header)).alias(id)])
                    .collect()?
                    .drop(&weights_header)?
            }
            let col_series = Series::new("Sequence_Name", &ids);
            dist_df.hstack(&[col_series])?
        }
        None => {
            let col_series = Series::new("Sequence_Name", &ids);
            dist_df.hstack(&[col_series])?
        }
    };

    Ok(dist_df)
}

pub fn compute_distance_matrix(
    fasta: &str,
    cluster_table: Option<&str>,
    label: &str,
    stringency: &Stringency,
    distance_method: &DistanceMethods,
) -> Result<()> {
    let (ids, sequences) = collect_fa_data(fasta)?;

    // double check that there are as many ids as there are sequences
    assert!(
        ids.len() == sequences.len(),
        "Unable to identify an ID for each sequence from the FASTA {}.",
        &fasta
    );

    // call a distance matrix with the chosen distance metric (defaulting to Levenshtein)
    let mut pw_distmat = SquareMatrix::from_pw_distances_with(&sequences, |seq1, seq2| {
        distance_method.calculate_distance(seq1, seq2)
    });
    pw_distmat.set_labels(ids.clone());

    // pull distmat information out of the SquareMatrix struct and convert to dataframe
    let mut dist_col_vec: Vec<Series> = vec![Default::default(); pw_distmat.size()];
    for (i, (column, label)) in pw_distmat
        .iter_cols()
        .zip(pw_distmat.iter_labels())
        .enumerate()
    {
        let series = Series::new(label, column.collect::<Vec<f64>>());
        dist_col_vec[i] = series;
    }

    let mut dist_df = process_cluster_info(cluster_table, dist_col_vec, &ids, stringency)?;

    // write out the weighted distance matrix
    let out_name = format!("{}-dist-matrix.csv", label);
    let out_handle = File::create(out_name).expect(
        "File could not be created to write the distance matrix to. Please check file-write permissions."
    );
    CsvWriter::new(out_handle)
        .has_header(true)
        .finish(&mut dist_df)
        .expect("Weighted distance matrix could not be written.");

    Ok(())
}
