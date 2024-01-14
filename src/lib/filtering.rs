use bio::io::fasta::Reader;
use chrono::NaiveDate;
use derive_new::new;
use noodles::{bgzf, fasta};
use polars::prelude::*;
use polars_io::ipc::IpcReader;
use std::collections::HashMap;
use std::fs::File;
use std::io::{self, BufRead, BufReader, Write};
use std::path::Path;

pub fn replace_gaps(input_path: &String) -> io::Result<()> {
    // define the input path
    let path: &Path = Path::new(input_path);

    // define the output file
    let path_out: &Path = Path::new("no-gaps.fasta");

    let file: File = File::open(path)?;
    let reader: BufReader<File> = BufReader::new(file);

    let mut output_file = File::create(path_out)?;

    // iterate through the FASTA records from the reader,
    // replace any "-" symbols with "N" characters, and
    // send them to the writer
    for line in reader.lines() {
        let mut line = line?;
        if !line.starts_with('>') {
            line = line.replace('-', "N");
        }
        writeln!(output_file, "{}", line)?;
    }

    Ok(())
}

#[derive(new)]
struct RefPriors {
    total_length: i32,
    max_n_count: f32,
}

fn determine_max_n(reference: &str, ambiguity: &f32) -> Result<RefPriors, std::io::Error> {
    let ref_handle = File::open(reference)?;
    let reader = Reader::new(ref_handle);

    let mut ref_length = 0;

    for record in reader.records() {
        ref_length += record.unwrap().seq().len() as i32;
    }

    let max_n_count = (ambiguity * ref_length as f32).floor();

    Ok(RefPriors::new(ref_length, max_n_count))
}

pub fn filter_by_n(input_path: &str, ambiguity: &f32, reference: &str) -> io::Result<()> {
    let mut fasta_reader = File::open(input_path)
        .map(bgzf::Reader::new)
        .map(fasta::Reader::new)?;

    // determine the maximum number of masked bases given the provided reference
    let ref_priors = match determine_max_n(reference, ambiguity) {
        Ok(value) => value,
        Err(message) => return Err(message),
    };

    let path_out: &Path = Path::new("filtered-by-n.fasta");
    let output_handle = File::create(path_out)?;
    let mut fasta_writer: fasta::Writer<File> = fasta::Writer::new(output_handle);

    for record in fasta_reader.records() {
        let record = record.expect("Error reading record");
        let sequence = record.sequence();
        let sequence_length = sequence.len();
        let id = record.name();

        // compute necessary information given the reference
        let length_diff = ref_priors.total_length - (sequence_length as i32);
        let incompleteness = if length_diff < 0 { 0 } else { length_diff };

        // Count the number of "N" characters in the sequence
        let count_n = (sequence
            .get(..)
            .unwrap()
            .iter()
            .filter(|&&c| c as char == 'N')
            .count() as i32
            + incompleteness) as f32;

        // write to output FASTA if N-count is less than the desired amount
        if count_n <= ref_priors.max_n_count {
            println!(
                "Record {} had more than {} masked bases and is thus removed.",
                id, ref_priors.max_n_count
            )
        } else {
            fasta_writer
                .write_record(&record)
                .expect("Error writing record");
        }
    }

    Ok(())
}

pub fn date_accessions(input_path: &str) -> Result<HashMap<String, String>, PolarsError> {
    // Read the Arrow file into a DataFrame
    let file = File::open(input_path).expect("Cannot open file");
    let df = IpcReader::new(file).finish()?;

    // Make sure the required columns exist
    if !df.schema().contains("Accession") || !df.schema().contains("Isolate Collection date") {
        return Err(PolarsError::ColumnNotFound(
            "Required columns not found".into(),
        ));
    }

    // Create a HashMap from the "Accession" and "Isolate Collection date" columns
    let accession_s = df
        .column("Accession")?
        .cast(&DataType::Utf8)
        .expect("Cannot cast to Utf8");
    let date_s = df
        .column("Isolate Collection date")?
        .cast(&DataType::Utf8)
        .expect("Cannot cast to Utf8");

    let mut accession_to_date = HashMap::new();
    for (accession, date) in accession_s
        .utf8()
        .unwrap()
        .into_iter()
        .zip(date_s.utf8().unwrap().into_iter())
    {
        if let (Some(accession), Some(date)) = (accession, date) {
            accession_to_date.insert(accession.to_string(), date.to_string());
        }
    }

    Ok(accession_to_date)
}

fn lookup_date(record_name: &str, lookup: &HashMap<String, String>) -> Option<NaiveDate> {
    match lookup.get(record_name) {
        Some(date) => match NaiveDate::parse_from_str(date, "%Y-%m-%d") {
            Ok(parsed_date) => Some(parsed_date),
            Err(_) => None,
        },
        None => None,
    }
}

pub fn separate_by_month(
    input_fasta: &str,
    accession_to_date: HashMap<String, String>,
) -> io::Result<()> {
    let mut fasta_reader = File::open(input_fasta)
        .map(bgzf::Reader::new)
        .map(fasta::Reader::new)?;

    let mut open_writers: HashMap<String, fasta::Writer<File>> = HashMap::new();

    for record in fasta_reader.records() {
        let record = record?;
        let record_date = lookup_date(record.name(), &accession_to_date);

        match record_date {
            Some(date) => {
                let year_month = date.format("%Y-%m").to_string();
                let writer = open_writers.entry(year_month.clone()).or_insert_with(|| {
                    let file = File::create(format!("{}.fasta", year_month)).unwrap();
                    fasta::Writer::new(file)
                });

                writer.write_record(&record)?;
            }
            None => println!(
                "Skipping record with missing or unparseable date: {}",
                record.name()
            ),
        }
    }

    Ok(())
}
