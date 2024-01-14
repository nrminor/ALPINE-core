pub mod lib {
    // use std::{fs::File, io::BufReader};

    // use noodles::{bgzf, fasta};

    pub mod distmat;
    pub mod filtering;
    pub mod prevalence;
    pub mod reporting;

    // pub enum SupportedCodec {
    //     Gzip,
    //     Zstd,
    //     Bzip2,
    //     BGzip,
    // }

    // trait DecodeInput {
    //     fn get_fasta_reader(
    //         &self,
    //         fasta_path: &str,
    //     ) -> anyhow::Result<fasta::Reader<fasta::Reader<std::fs::File>>>;
    // }

    // impl DecodeInput for SupportedCodec {
    //     fn get_fasta_reader(
    //         &self,
    //         fasta_path: &str,
    //     ) -> anyhow::Result<fasta::Reader<fasta::Reader<std::fs::File>>> {
    //         match &self {
    //             SupportedCodec::Gzip => (),
    //             SupportedCodec::Zstd => (),
    //             SupportedCodec::Bzip2 => (),
    //             SupportedCodec::BGzip => {
    //                 let mut fasta_reader = File::open(fasta_path)
    //                     .map(bgzf::Reader::new)
    //                     .map(fasta::Reader::new);
    //                 fasta_reader
    //             }
    //             _ => {
    //                 let mut fasta_reader = File::open(fasta_path)
    //                     .map(BufReader::new)
    //                     .map(fasta::Reader::new);
    //                 fasta_reader
    //             }
    //         }
    //     }
    // }
}
