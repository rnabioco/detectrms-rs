use rust_htslib::bam::ext::BamRecordExtensions;
use rust_htslib::bam::pileup::Indel;
use rust_htslib::bam::Read;
use seq_io::fasta::{Reader as FastaReader, Record};
use std::fs::File;
use std::io::{BufReader, Write, BufWriter};
use std::io::{self};
use std::env;
use std::process;
use std::collections::HashMap;
use rayon::prelude::*;
use std::sync::{Arc, Mutex};
use std::process::Command;
use std::str;
use rust_htslib::bam::IndexedReader;
use rust_htslib::bam::FetchDefinition;

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
    let mut bam_reader = IndexedReader::from_path(bam_file).expect("Error opening BAM file");
    let mut fasta_reader = FastaReader::new(BufReader::new(File::open(fasta_file).expect("Error opening FASTA file")));

    let tid = bam_reader.header().tid(ref_name.as_bytes()).expect("Reference not found") as i32;
    let start: i64 = 0;
    let end: i64 = 1 << 29; // Large value to cover the range of a reference

    bam_reader.fetch(FetchDefinition::Region(tid, start, end)).expect("Error setting region");

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
    //println!("{}: {} reads, coverage {}-{}", ref_name, ref_reads, ref_cov_start, ref_cov_end);
    //prints  8396 reads, coverage 0-18761 for every reference
    if ref_reads < 10 {
        println!("Skipping reference {} with only {} reads", ref_name, ref_reads);
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
        println!("Processing position {}", pos + 1);  // Debug print, this is called


        // Fetch pileup for the current position
        for pileup in bam_reader.pileup() {
                  
            println!("Processing pileup at positionxs");//this is never called
            let pileup_column = pileup.expect("Error reading pileup column");
            println!("{:?}", pileup_column.pos());
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
        println!("{}:{}", pos + 1, coverage);  // Debug print
    
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

    // Run `samtools idxstats` and capture the output
    let output = Command::new("samtools")
        .arg("idxstats")
        .arg(bam_file)
        .output()
        .expect("Failed to execute samtools");

    // Convert the output to a string
    let output_str = str::from_utf8(&output.stdout).expect("Failed to read stdout");

    // Parse the output to get references with reads
    let references_with_reads: Vec<String> = output_str.lines()
        .filter_map(|line| {
            let fields: Vec<&str> = line.split('\t').collect();
            if fields.len() >= 3 && fields[2] != "0" {
                Some(fields[0].to_string())
            } else {
                None
            }
        })
        .collect();

    println!("Number of bam references with reads: {}", references_with_reads.len());

    // Create and set up the output file with Arc and Mutex for thread safety
    let output_file = File::create(output_file).expect("Error creating output file");
    let output_file = Arc::new(Mutex::new(BufWriter::new(output_file)));

    // Write header to output file
    {
        let mut output = output_file.lock().unwrap();
        writeln!(output, "#Ref\tpos\tbase\tstrand\tcov\tmean_q\tmedian_q\tstd_q\tmis\tins\tdel\tACGT").expect("Error writing header to output file");
    }

    // Parallel processing of filtered references
    references_with_reads.par_iter().for_each(|ref_name_str| {
        println!("Processing reference {}...", ref_name_str);
        let results = process_ref(bam_file, fasta_file, ref_name_str);

        // Concurrently write results to the file
        let mut output = output_file.lock().unwrap();
        for line in results {
            println!("Writing line: {}", line);  // Debug print
            writeln!(output, "{}", line).expect("Error writing to output file");
        }
    });
}