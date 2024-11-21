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
    fastqualitydrop(
        &args.reads_1_arg,
        &args.reads_2_arg,
    );

     fastq_quality_drop(
        &args.reads_1_arg,
        &args.reads_2_arg,
    );
}

fn qualityscore() -> (Vec<usize>, Vec<&'static str>) {
    let qualitystring = (33..75).collect::<Vec<usize>>();
    let qualitydrop: Vec<_> = vec!["!", "\"\"", "#", "$", "%", "&","\'", "(", ")", "*", "+",",", "-", ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":", ";", "<", "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"];
    (qualitystring, qualitydrop)
}

fn fastqualitydrop(fastq1: &str, fastq2: &str){

     let (qualitystring, qualitydrop) = qualityscore();

}

fn fastq_quality_drop(fastq1: &str, fastq2: &str){

    let (qualitystring, qualitydrop) = qualityscore();
}
