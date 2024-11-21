mod args;
use args::FastqArgs;
use clap::Parser;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

/*
 *Author Gaurav Sablok
 *Universitat Potsdam
 *Date 2024-11-22

* implementing all the parts of the fastp in rustlang for the faster read user
* in nextseqseq, novaseq, and other high-throughput illumina platforms. These are
* available as separate rust-applications each of them and also a combined rust-fastp.
* This is the rust-fastp-quality-drop, which will drop the read bases at the specific position.
*
* */

fn main() {
    let args = FastqArgs::parse();
    fastqqualitydrop(
        &args.reads_1_arg,
        &args.reads_2_arg,
    );

     fastq_quality_drop(
        &args.reads_1_arg,
        &args.reads_2_arg,
    );
}

