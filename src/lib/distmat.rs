use anyhow::anyhow;
use anyhow::Result;
use clap::ValueEnum;
use derive_new::new;
use distmat::SquareMatrix;
use noodles::{bgzf, fasta};
use polars::{lazy::dsl::col, prelude::*};
use std::fmt;
use std::fs::File;
use std::io::BufReader;
use std::io::ErrorKind;
use std::ops::Mul;
use textdistance::{
    nstr::{lcsseq, lcsstr},
    str::{damerau_levenshtein, jaro_winkler, levenshtein, ratcliff_obershelp, smith_waterman},
    str::{hamming, jaccard},
};

#[derive(ValueEnum, Debug, Clone, PartialEq)]
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

impl fmt::Display for DistanceMethods {
    fn fmt(&self, f: &mut fmt::Formatter<'_>) -> fmt::Result {
        write!(
            f,
            "{}",
            match self {
                DistanceMethods::Hamming => "hamming",
                DistanceMethods::Levenshtein => "levenshtein",
                DistanceMethods::DamerauLevenshtein => "damerau-levenshtein",
                DistanceMethods::JaroWinkler => "jaro-winkler",
                DistanceMethods::SmithWaterman => "smith-waterman",
                DistanceMethods::RatcliffObershelp => "ratcliff-obershelp",
                DistanceMethods::LCSSeq => "lcs-seq",
                DistanceMethods::LCSStr => "lcs-str",
                DistanceMethods::Jaccard => "jaccard",
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

/// Cluster columns contains the column names where the information ALPINE needs is stored
#[derive(new, Debug, Clone)]
struct ClusterColumns {
    type_col: String,
    index_col: String,
    id_col: String,
    size_col: String,
}

fn get_cluster_cols(cluster_table: &LazyFrame) -> Result<ClusterColumns> {
    // separate out the columns of information we need
    let cluster_query = cluster_table.clone().collect()?;
    let col_names = cluster_query.get_column_names();

    let type_col = match col_names.first() {
        Some(col_name) => col_name.to_string(),
        None => {
            eprintln!(
                "Please double check that the column of VSEARCH cluster types is the first column."
            );
            return Err(anyhow!(
                "Member types could not be parsed from provided cluster table,"
            ));
        }
    };

    let index_col = match col_names.get(1) {
        Some(col_name) => col_name.to_string(),
        None => {
            eprintln!(
                "Please double check that the column of VSEARCH cluster index is the second column."
            );
            return Err(anyhow!(
                "Column indices could not be parsed from provided cluster table,"
            ));
        }
    };

    let name_col = match col_names.get(8) {
        Some(col_name) => col_name.to_string(),
        None => {
            eprintln!("Please double check that the column of sequence names is the ninth column.");
            return Err(anyhow!(
                "Sequence names could not be parsed from provided cluster table,"
            ));
        }
    };

    let size_col = match col_names.get(2) {
        Some(col_name) => col_name.to_string(),
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

fn compute_weight(
    centroid_lf: LazyFrame,
    clust_index: i64,
    index_col: &str,
    size_col: &str,
) -> Result<f64> {
    let collected_df = centroid_lf
        .clone()
        .filter(col(index_col).eq(clust_index))
        .select([col(size_col)])
        .collect()?;

    let attempt = match collected_df
        .column(size_col)?
        .iter()
        .next() {
            Some(value) => value,
            None => return Err(anyhow!("Could not parse centroid data to compute a weight. Please double check the input cluster table."))
        };

    Ok(attempt.try_extract::<f64>()?)
}

fn weight_by_cluster_size(
    seq_name: &str,
    stringency: &Stringency,
    cluster_table: &LazyFrame,
) -> Result<(String, Series)> {
    let clust_cols = get_cluster_cols(cluster_table)?;

    // Filter down the df so that only rows representing centroids are present
    let centroid_lf = cluster_table
        .clone()
        .filter(col(&clust_cols.type_col).eq(lit("C")))
        .sort(&clust_cols.id_col, Default::default());

    // Filter down to hits only and use to get a total number of sequences for the current month
    let hits_df = cluster_table
        .clone()
        .filter(
            col(&clust_cols.type_col)
                .eq(lit("H"))
                .or(col(&clust_cols.type_col).eq(lit("S"))),
        )
        .sort(&clust_cols.id_col, Default::default())
        .collect()?;
    let month_total: f64 = if hits_df.shape().0 == 0 {
        1.0
    } else {
        hits_df.shape().0 as f64
    };

    // determine the cluster index for the current cluster member
    let index: i64 = cluster_table
        .clone()
        .filter(col(&clust_cols.id_col).eq(lit(seq_name)))
        .select(&[col(&clust_cols.index_col)])
        .collect()?
        .column(&clust_cols.index_col)?
        .get(0)?
        .try_extract::<i64>()?;

    // find cluster size for the member accession in question
    let weighting_size = compute_weight(
        centroid_lf,
        index,
        &clust_cols.index_col,
        &clust_cols.size_col,
    )?;
    let weighting_freq = weighting_size / month_total;

    // compute weights lazyframe column
    let weights_header = format!("{}_weights", seq_name);
    let weights_lf = match *stringency {
        Stringency::Strict | Stringency::Extreme => hits_df
            .lazy()
            .with_column(
                ((col(&clust_cols.size_col) * lit(weighting_freq.ln() * -1.0)) / lit(month_total))
                    .alias(&weights_header),
            )
            .select(&[col(&weights_header)]),
        _ => hits_df
            .lazy()
            .with_column(((lit(1) - col(&clust_cols.size_col)) / lit(month_total)).alias("Weights"))
            .select(&[col(&weights_header)]),
    };

    let weights = weights_lf.collect()?.column(&weights_header)?.to_owned();

    // double check that we now have as many weights as we need
    assert!(
        weights.len() as f64 == month_total,
        "ALPINE was not able to find the same number of clusters as are represented in the centroid distance matrix"
    );

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

fn collect_fa_data(fasta: &str) -> Result<(Vec<String>, Vec<String>)> {
    // pull out record IDs and sequences into their own string vectors, while
    // handling potential bgzip compression
    let parsed_fasta: std::io::Result<Vec<(String, String)>> = if fasta.ends_with(".gz") {
        File::open(fasta)
            .map(bgzf::Reader::new)
            .map(fasta::Reader::new)?
            .records()
            .map(|result| {
                result.and_then(|record| {
                    let id = record.name().to_owned();
                    unpack_sequence(&record).map(|sequence_string| (id, sequence_string))
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
                    unpack_sequence(&record).map(|sequence_string| (id, sequence_string))
                })
            })
            .collect()
    };

    let (ids, sequences) = match parsed_fasta {
        Ok(pairs) => pairs.into_iter().unzip(),
        Err(e) => return Err(e.into()),
    };

    Ok((ids, sequences))
}

pub fn compute_distance_matrix(
    fasta: &str,
    cluster_table: &str,
    yearmonth: &str,
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
    let mut dist_col_vev: Vec<Series> = vec![Default::default(); pw_distmat.size()];
    for (i, (column, label)) in pw_distmat
        .iter_cols()
        .zip(pw_distmat.iter_labels())
        .enumerate()
    {
        let series = Series::new(label, column.collect::<Vec<f64>>());
        dist_col_vev[i] = series;
    }
    let mut dist_df = DataFrame::new(dist_col_vev)?;

    // read the cluster table into a lazyframe to query for cluster-size-based distance weights
    let cluster_df = CsvReader::from_path(cluster_table)?
        .has_header(false)
        .with_delimiter(b'\t')
        .finish()?
        .lazy();

    // multiply weights onto each column based on the sequence it represents
    for id in ids {
        let (weights_header, weights) = weight_by_cluster_size(&id, stringency, &cluster_df)?;
        dist_df = dist_df
            .hstack(&[weights])?
            .lazy()
            .with_columns(&[col(&id).mul(col(&weights_header)).alias(&id)])
            .collect()?
            .drop(&weights_header)?
    }

    // write out the weighted distance matrix
    let out_name = format!("{}-distmat.csv", yearmonth);
    let out_handle = File::create(out_name).expect(
        "File could not be created to write the distance matrix to. Please check file-write permissions."
    );
    CsvWriter::new(out_handle)
        .has_header(true)
        .finish(&mut dist_df)
        .expect("Weighted distance matrix could not be written.");

    Ok(())
}
