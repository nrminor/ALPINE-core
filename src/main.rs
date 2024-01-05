use alpine::lib::distmat::{DistanceMethods, Stringency};
use alpine::lib::*;
use anyhow::Result;
use clap::{Parser, Subcommand};
use tokio::runtime::Builder;

const INFO: &str = "ALPINE Rust command line interface.";

#[derive(Parser)]
#[clap(name = "alpine")]
#[clap(about = INFO)]
struct Cli {
    #[arg(short, long, action = clap::ArgAction::Count)]
    debug: u8,

    #[command(subcommand)]
    command: Option<Commands>,
}

#[derive(Subcommand)]
enum Commands {
    ReplaceGaps {
        #[arg(required = true)]
        fasta: String,
    },

    FilterByNs {
        #[arg(required = true)]
        fasta: String,
        ambiguity: f32,
    },

    SeparateByMonth {
        #[arg(required = true)]
        fasta: String,
        metadata: String,
    },

    DistanceMatrix {
        #[arg(required = true)]
        fasta: String,
        cluster_table: String,
        yearmonth: String,
        stringency: Stringency,
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
        Some(Commands::ReplaceGaps { fasta }) => {
            filtering::replace_gaps(fasta)?;
            Ok(())
        }
        Some(Commands::FilterByNs { fasta, ambiguity }) => {
            filtering::filter_by_n(fasta, ambiguity)?;
            Ok(())
        }
        Some(Commands::SeparateByMonth { fasta, metadata }) => {
            let accession_to_date = filtering::date_accessions(metadata)?;
            filtering::separate_by_month(fasta, accession_to_date)?;
            Ok(())
        }
        Some(Commands::DistanceMatrix { .. }) => {
            distmat::compute_distance_matrix()?;
            Ok(())
        }
        None => {
            println!("{}\n", INFO);
            std::process::exit(1);
        }
    }
}

fn main() {
    let ncores = 4;

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
