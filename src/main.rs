use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::{Reader as BamReader, Read};
use seq_io::fasta::{Reader as FastaReader, Record};
use std::fs::File;
use std::io::{BufReader, Write, BufWriter};
use std::io::{self};
use std::env;
use std::process;
use std::collections::HashMap;

fn fetch_reference_base(fasta_reader: &mut FastaReader<BufReader<File>>, ref_name: &str, pos: usize) -> io::Result<char> {
    for result in fasta_reader.records() {
        let record = match result {
            Ok(record) => record,
            Err(e) => return Err(io::Error::new(io::ErrorKind::Other, format!("Fasta read error: {}", e))),
        };
        if let Ok(record_id) = record.id() {
            if record_id == ref_name {
                let sequence = record.seq();
                if pos < sequence.len() {
                    return Ok(sequence[pos] as char);
                } else {
                    break; // Position is out of range for this reference
                }
            }
        }
    }
    Err(io::Error::new(io::ErrorKind::NotFound, "Reference or position not found"))
}  

fn process_ref(bam_file: &str, fasta_file: &str, ref_name: &str) -> Vec<String> {
    let mut bam_reader = BamReader::from_path(bam_file).expect("Error opening BAM file");
    let mut fasta_reader = FastaReader::new(BufReader::new(File::open(fasta_file).expect("Error opening FASTA file")));

    let mut ref_reads = 0;
    let mut ref_cov_start = usize::MAX;
    let mut ref_cov_end = 0;
    let mut output_lines = Vec::new();

    for read_result in bam_reader.records() {
        let read = read_result.expect("Error reading BAM record");
        ref_reads += 1;
        //println!("{:?}", &ref_reads);
    
        let first_position = read.reference_start();
        let last_position = read.reference_end();
    
        if let Ok(first_pos_usize) = first_position.try_into() {
            if first_pos_usize < ref_cov_start {
                ref_cov_start = first_pos_usize;
            }
        }
    
        if let Ok(last_pos_usize) = last_position.try_into() {
            if last_pos_usize > ref_cov_end {
                ref_cov_end = last_pos_usize;
            }
        }
    }
    
    if ref_reads < 10 {
        eprintln!("Skipping reference {} with only {} reads", ref_name, ref_reads);
        return Vec::new();
    }

    // Process each position
    for pos in ref_cov_start..ref_cov_end {
        let mut base_counts = HashMap::from([('A', 0), ('C', 0), ('G', 0), ('T', 0)]);
        let mut qualities = Vec::new();
        let mut mismatches = 0;
        let mut insertions = 0;
        let mut deletions = 0;
        let mut strand = '+';
        let mut coverage = 0;

        // Fetch pileup for the current position
        for pileup in bam_reader.pileup() {
            let pileup_column = pileup.expect("Error reading pileup column");
    
            if pileup_column.pos() as usize != pos {
                continue;
            }
    
            for pileup_read in pileup_column.alignments() {
                coverage += 1;
                let is_del = pileup_read.is_del() == true;
                let is_refskip = pileup_read.is_refskip() == true;
                
                // Skip deletions and reference skips
                if is_del || is_refskip {
                    if is_del {
                        deletions += 1;
                    }
                    continue;
                }
    
                // Handle insertions
                match pileup_read.indel() {
                    Indel::Ins(_) => insertions += 1,
                    _ => {},
                }
    
                if let Some(qpos) = pileup_read.qpos() {
                    let base = pileup_read.record().seq()[qpos] as char;
                    let quality = pileup_read.record().qual()[qpos];
    
                    // Update base counts
                    *base_counts.entry(base).or_insert(0) += 1;
    
                    // Collect quality
                    qualities.push(quality);
    
                    // Determine strand
                    strand = if pileup_read.record().is_reverse() { '-' } else { '+' };
    
                    // Fetch the reference base and check for mismatches
                    // Assuming you have implemented a function to fetch the reference base
                    let reference_base = fetch_reference_base(&mut fasta_reader, ref_name, pos);
                    if base != reference_base.unwrap() {
                        mismatches += 1;
                    }
                }
            }
    

        // Calculate coverage
        coverage = base_counts.values().sum::<i32>();

        // Calculate mean quality
        let mean_quality = if !qualities.is_empty() {
            qualities.iter().sum::<u8>() as f64 / qualities.len() as f64
        } else {
            0.0
        };

        // Calculate median quality
        let median_quality = if !qualities.is_empty() {
            let mut sorted_qualities = qualities.clone();
            sorted_qualities.sort_unstable();
            let mid = sorted_qualities.len() / 2;
            if sorted_qualities.len() % 2 == 0 {
                (sorted_qualities[mid - 1] as f64 + sorted_qualities[mid] as f64) / 2.0
            } else {
                sorted_qualities[mid] as f64
            }
        } else {
            0.0
        };

        // Calculate standard deviation of quality
        let std_quality = if !qualities.is_empty() {
            let mean = mean_quality;
            let variance = qualities.iter().map(|&q| {
                let diff = q as f64 - mean;
                diff * diff
            }).sum::<f64>() / qualities.len() as f64;
            variance.sqrt()
        } else {
            0.0
        };

        // Calculate ratios
        let mismatch_ratio = if coverage > 0 { mismatches as f64 / coverage as f64 } else { 0.0 };
        let insertion_ratio = if coverage > 0 { insertions as f64 / coverage as f64 } else { 0.0 };
        let deletion_ratio = if coverage > 0 { deletions as f64 / coverage as f64 } else { 0.0 };

        // Format and append the output line
        output_lines.push(format!("{}\t{}\t{}\t{}\t{}\t{:.?}\t{:.?}\t{:.?}\t{:.?}\t{:.?}\t{:.?}\t{}:{}:{}:{}\n",
            ref_name, pos + 1, "reference_base", strand, coverage, mean_quality, median_quality, std_quality, mismatch_ratio, insertion_ratio, deletion_ratio,
            base_counts[&'A'], base_counts[&'C'], base_counts[&'G'], base_counts[&'T']));
    }
    }
    output_lines
}



fn main() {
    let args: Vec<String> = env::args().collect();

    if args.len() != 4 {
        eprintln!("Usage: {} <bam_file> <fasta_file> <output_file>", args[0]);
        process::exit(1);
    }

    let bam_file = &args[1];
    let fasta_file = &args[2];
    let output_file = &args[3];

    let bam_reader = BamReader::from_path(bam_file).expect("Error opening BAM file");

    println!("Creating output file at: {}", output_file);

    let mut output_file = BufWriter::new(File::create(output_file).expect("Error creating output file"));

    println!("Writing header to output file");

    if let Err(e) = writeln!(output_file, "#Ref\tpos\tbase\tstrand\tcov\tmean_q\tmedian_q\tstd_q\tmis\tins\tdel\tACGT") {
        eprintln!("Error writing header to output file: {}", e);
        return;
    }

    output_file.flush().expect("Error flushing output file");

    // Process each reference
    for ref_name_bytes in bam_reader.header().target_names() {
        let ref_name_str = std::str::from_utf8(ref_name_bytes).expect("Found invalid UTF-8");
        println!("Processing reference {}...", ref_name_str);
        let results = process_ref(bam_file, fasta_file, ref_name_str);

        // Write results to the file
        for line in results {
            writeln!(output_file, "{}", line).expect("Error writing to output file");
        }
        output_file.flush().expect("Error flushing output file");

    }
}