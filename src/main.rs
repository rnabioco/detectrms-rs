use clap::Parser;
use rust_htslib::bam::{Read, IndexedReader};
use std::collections::{HashMap, BTreeMap};
use std::path::Path;
use bio::io::fasta;
use std::fs::File;
use std::io::{Write, BufWriter};

fn mean(qualities: &[u8]) -> f64 {
    if qualities.is_empty() {
        0.0
    } else {
        qualities.iter().map(|&q| q as f64).sum::<f64>() / qualities.len() as f64
    }
}

fn median(qualities: &[u8]) -> f64 {
    if qualities.is_empty() {
        return 0.0;
    }

    let mut sorted = qualities.to_vec();
    sorted.sort_unstable();
    let mid = sorted.len() / 2;
    
    if sorted.len() % 2 == 0 {
        (sorted[mid - 1] as f64 + sorted[mid] as f64) / 2.0
    } else {
        sorted[mid] as f64
    }
}

fn std_dev(qualities: &[u8], mean: f64) -> f64 {
    let variance = qualities
        .iter()
        .map(|&q| {
            let diff = q as f64 - mean;
            diff * diff
        })
        .sum::<f64>() / qualities.len() as f64;
    variance.sqrt()
}

#[derive(Parser, Debug)]
#[clap(author, version, about, long_about = None)]
struct Args {
    /// Input BAM file
    #[clap(short, long, value_parser)]
    bam: String,

    /// Input FASTA file
    #[clap(short, long, value_parser)]
    fasta: String,

    /// Output TSV file
    #[clap(short, long, value_parser)]
    output: String,
}

fn main() {
    let args = Args::parse();

    let bam_path = Path::new(&args.bam);
    let mut bam_reader = IndexedReader::from_path(&bam_path).expect("Error opening BAM file");
    let fasta_path = Path::new(&args.fasta);
    let fasta_reader = fasta::Reader::from_file(&fasta_path).expect("Error opening FASTA file");

    // Prepare to write to a TSV file
    let output_file = File::create(&args.output).expect("Error creating output file");
    let mut output = BufWriter::new(output_file);
    // Write the header
    writeln!(output, "#Ref\tpos\tbase\tstrand\tcov\tmean_q\tmedian_q\tstd_q\tmis\tins\tdel\tACGT").expect("Error writing header");

    // Load sequences from FASTA
    let mut sequences = BTreeMap::new();
    for result in fasta_reader.records() {
        let record = result.expect("Error during fasta record parsing");
        let tid = record.id().to_string();
        sequences.insert(tid, record.seq().to_owned());
    }

    // Iterate over sequences and process pileups
    for (tid, sequence) in &sequences {
        match bam_reader.fetch((tid, 0, 1 << 29)) {
            Ok(_) => {},
            Err(e) => eprintln!("Error setting region for TID: {}: {}", tid, e),
        }

        for pileup in bam_reader.pileup() {
            if let Ok(p) = pileup {
                let pos = p.pos() as usize;
                let depth = p.depth();
                let mut base_counts = HashMap::from([('A', 0), ('C', 0), ('G', 0), ('T', 0)]);
                let mut qualities: Vec<u8> = Vec::new();
                let mut mismatches = 0;
                let mut insertions = 0;
                let mut deletions = 0;
                let mut strand: char = ' ';
                let reference_base = *sequence.get(pos).unwrap_or(&b'N') as char;

                for pileup_read in p.alignments() {
                    if pileup_read.record().is_reverse() {
                        strand = '-';
                    } else {
                        strand = '+';
                    }
                    let is_del = pileup_read.is_del();
                    let is_refskip = pileup_read.is_refskip();

                    if is_del || is_refskip {
                        if is_del {
                            deletions += 1;
                        }
                        continue;
                    }

                    if let Some(qpos) = pileup_read.qpos() {
                        let base = pileup_read.record().seq()[qpos] as char;
                        let quality = pileup_read.record().qual()[qpos];

                        *base_counts.entry(base).or_insert(0) += 1;
                        qualities.push(quality);

                        if base != reference_base {
                            mismatches += 1;
                        }
                    }
                }

                let coverage = base_counts.values().sum::<i32>();
                let mean_quality = mean(&qualities);
                let median_quality = median(&qualities);
                let std_quality = std_dev(&qualities, mean_quality);
                let mismatch_ratio = mismatches as f64 / depth as f64;
                let insertion_ratio = insertions as f64 / depth as f64;
                let deletion_ratio = deletions as f64 / depth as f64;

                let line = format!(
                    "{}\t{}\t{}\t{}\t{}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{:.2}\t{}:{}:{}:{}",
                    tid,
                    pos + 1,
                    reference_base,
                    strand,
                    coverage,
                    mean_quality,
                    median_quality,
                    std_quality,
                    mismatch_ratio,
                    insertion_ratio,
                    deletion_ratio,
                    base_counts[&'A'],
                    base_counts[&'C'],
                    base_counts[&'G'],
                    base_counts[&'T']
                );

                writeln!(output, "{}", line).expect("Error writing to file");
            } else {
                eprintln!("Error reading pileup");
            }
        }
    }
}