use alpine::lib::distmat::{DistanceMethods, Stringency};
use alpine::lib::*;
use anyhow::Result;
use clap::{Parser, Subcommand};
use tokio::runtime::Builder;

/// Shout out to https://patorjk.com/software/taag/ for the ASCII art.
const INFO: &str = r"

_____/\\\\\\\\\_____/\\\______________/\\\\\\\\\\\\\____/\\\\\\\\\\\__/\\\\\_____/\\\__/\\\\\\\\\\\\\\\_
 ___/\\\\\\\\\\\\\__\/\\\_____________\/\\\/////////\\\_\/////\\\///__\/\\\\\\___\/\\\_\/\\\///////////__
  __/\\\/////////\\\_\/\\\_____________\/\\\_______\/\\\_____\/\\\_____\/\\\/\\\__\/\\\_\/\\\_____________
   _\/\\\_______\/\\\_\/\\\_____________\/\\\\\\\\\\\\\/______\/\\\_____\/\\\//\\\_\/\\\_\/\\\\\\\\\\\_____
    _\/\\\\\\\\\\\\\\\_\/\\\_____________\/\\\/////////________\/\\\_____\/\\\\//\\\\/\\\_\/\\\///////______
     _\/\\\/////////\\\_\/\\\_____________\/\\\_________________\/\\\_____\/\\\_\//\\\/\\\_\/\\\_____________
      _\/\\\_______\/\\\_\/\\\_____________\/\\\_________________\/\\\_____\/\\\__\//\\\\\\_\/\\\_____________
       _\/\\\_______\/\\\_\/\\\\\\\\\\\\\\\_\/\\\______________/\\\\\\\\\\\_\/\\\___\//\\\\\_\/\\\\\\\\\\\\\\\_
        _\///________\///__\///////////////__\///______________\///////////__\///_____\/////__\///////////////__

ALPINE: Anachronistic Lineage and Persistent INfection Explorer
===============================================================

Command line interface for the core Rust components of ALPINE.
These commands are called in the full pipeline, which is imple-
mented in Nextflow alongside `seqkit`, `csvtk`, `nushell`,
`vsearch`, and a bin of bespoke Python scripts. However, users
may also use the commands available in this crate to run
similar analyses themselves via the command line.
";

#[derive(Parser)]
#[clap(name = "alpine")]
#[clap(about = INFO)]
#[clap(version = "v0.1.2")]
struct Cli {
    #[command(flatten)]
    verbose: clap_verbosity_flag::Verbosity,

    #[arg(short, long, default_value_t = 3, required = false)]
    threads: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    #[clap(
        about = "Replace FASTA gap symbols '-' with masked bases 'N'.",
        aliases = &["rg", "gap"]
    )]
    ReplaceGaps {
        /// FASTA format file of sequences. Reads from stdin if no path is provided.
        #[arg(short, long, required = false)]
        fasta: Option<String>,

        /// Output file path. Writes to stdout if no path is provided.
        #[arg(short, long, required = false)]
        output_file: Option<String>,
    },

    #[clap(
        about = "Filter out FASTA records with more than the desired number of masked 'N' bases.",
        aliases = &["fbn", "mask", "nbase"]
    )]
    FilterByN {
        /// FASTA format file of sequences
        #[arg(short, long, required = true)]
        fasta: String,

        // Maximum number of ambiguous masked bases ("N") to tolerate in each FASTA record.
        #[arg(short, long, required = true)]
        ambiguity: f32,

        /// Reference sequence in FASTA format for the taxon of interest.
        #[arg(short, long, required = true)]
        reference: String,

        /// Output file path. Writes to stdout if no path is provided.
        #[arg(short, long, required = false)]
        out_file: Option<String>,
    },

    #[clap(
        about = "Use collection dates from FASTA record metadata to sort all FASTA records into a separate FASTA for each year-month combination.",
        aliases = &["month", "sbm"]
    )]
    SeparateByMonth {
        /// FASTA format file of sequences. Reads from stdin if no path is provided.
        #[arg(short, long, required = false)]
        fasta: Option<String>,

        /// Output file path. Writes to stdout if no path is provided.
        #[arg(short, long, required = true)]
        metadata: String,
    },

    #[clap(
        about = "Compute a symmetric pairwise distance matrix based on how dissimilar sequences in the provided FASTA are to one another.",
        aliases = &["distmat", "dm"]
    )]
    DistanceMatrix {
        /// FASTA format file of sequences
        #[arg(short, long, required = true)]
        fasta: String,

        /// Tab-delimited table of cluster metadata provided by `VSEARCH` clustering algorithms.
        #[arg(short, long, required = false)]
        cluster_table: Option<String>,

        /// Label for the output distance matrix. The ALPINE pipeline uses a chronological "year-date" here.
        #[arg(short, long, required = true)]
        label: String,

        /// How strictly to differentiate high- from low-distance clusters.
        #[arg(short, long, required = false, default_value_t = Stringency::Strict)]
        stringency: Stringency,

        /// Choice of distance-calling methods.
        #[arg(short, long, default_value_t = DistanceMethods::BlockAlignment)]
        distance_method: DistanceMethods,
    },
}

pub fn print_errors(response: Result<()>) {
    match response {
        Ok(_) => {}
        Err(err) => eprintln!("Error: {}", err),
    }
}

async fn run() -> Result<()> {
    let cli = Cli::parse();
    match &cli.command {
        Some(Commands::ReplaceGaps { fasta, output_file }) => {
            filtering::replace_gaps(fasta.as_deref(), output_file.as_deref())?;
            Ok(())
        }
        Some(Commands::FilterByN {
            fasta,
            ambiguity,
            reference,
            out_file,
        }) => {
            filtering::filter_by_n(fasta, ambiguity, reference, out_file.as_deref())?;
            Ok(())
        }
        Some(Commands::SeparateByMonth { fasta, metadata }) => {
            let accession_to_date = filtering::date_accessions(metadata)?;
            filtering::separate_by_month(fasta.as_deref(), accession_to_date)?;
            Ok(())
        }
        Some(Commands::DistanceMatrix {
            fasta,
            cluster_table,
            label,
            stringency,
            distance_method,
        }) => {
            distmat::compute_distance_matrix(
                fasta,
                cluster_table.as_deref(),
                label,
                stringency,
                distance_method,
            )?;
            Ok(())
        }
        None => {
            eprintln!("{}\n", INFO);
            std::process::exit(1);
        }
    }
}

fn main() {
    // Parse the number of threads to use
    let ncores = Cli::parse().threads as usize;

    let runtime = Builder::new_multi_thread()
        .worker_threads(ncores)
        .enable_all()
        .build()
        .unwrap();

    runtime.block_on(async {
        match run().await {
            Ok(_) => {}
            Err(e) => {
                eprintln!("Error: {:?}", e);
                std::process::exit(1);
            }
        }
    });
}
