mod args;
use std::str::Chars;
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
    fastq_quality_drop(
        &args.reads_1_arg,
        &args.reads_2_arg,
        args.quality_score,
    );

}

fn qualityscore() -> (Vec<usize>, Vec<Chars<'static>>) {
   let dropval: Vec<_> = vec!["!", "\"", "#", "$", "%", "&","\'", "(", ")", "*", "+",
    ",", "-", ".", "/", "0", "1", "2", "3", "4", "5", "6", "7", "8", "9", ":", ";", "<",
    "=", ">", "?", "@", "A", "B", "C", "D", "E", "F", "G", "H", "I", "J", "K"]
    .iter()
    .map(|x| x.chars())
    .collect::<Vec<_>>();
   let dropvec = (0..43).collect::<Vec<usize>>();
    (dropvec, dropval)
}

fn fastq_quality_drop(fastq1: &str, fastq2: &str, threshold: usize){

     let (qualitystring, qualitydrop) = qualityscore();

    #[derive(Debug, Clone)]
    struct FileFastqPre {
     header: String,
     sequence: String,
     strand:String,
     quality:String,
    }


    let file1 = File::open(fastq1).expect("file not present");
    let file2 = File::open(fastq2).expect("file not present");
    let fileread_1 = BufReader::new(&file1);
    let fileread_2 = BufReader::new(&file2);
    let mut fastq_1: Vec<FileFastqPre> = Vec::new();
    let mut fastq_2: Vec<FileFastqPre> = Vec::new();
    let mut fastq_1_header: Vec<String> = Vec::new();
    let mut fastq_2_header: Vec<String> = Vec::new();
    let mut fastq_1_sequence: Vec<String> = Vec::new();
    let mut fastq_2_sequence: Vec<String> = Vec::new();
    let mut fastq_1_strand: Vec<String> = Vec::new();
    let mut fastq_2_strand: Vec<String> = Vec::new();
    let mut fastq_1_quality: Vec<String> = Vec::new();
    let mut fastq_2_quality: Vec<String> = Vec::new();


    for i in fileread_1.lines(){
    let line = i.expect("line not present");
    if line.starts_with("@"){
        fastq_1_header.push(line);
    } else if line.starts_with("A") && !line.contains("E") || line.starts_with("T")
        && !line.contains("E") || line.starts_with("G") && !line.contains("E") ||
            line.starts_with("C") && !line.contains("E") || line.starts_with("N") && !line.contains("E")
    {
    fastq_1_sequence.push(line);
    } else if line.starts_with("+") || line.starts_with("-") {
       fastq_1_strand.push(line);
    } else if line.contains("E"){
       fastq_1_quality.push(line);
    }
  }

    for i in fileread_2.lines(){
    let line = i.expect("line not present");
    if line.starts_with("@"){
        fastq_2_header.push(line);
    } else if line.starts_with("A") && !line.contains("E") || line.starts_with("T")
        && !line.contains("E") || line.starts_with("G") && !line.contains("E") ||
            line.starts_with("C") && !line.contains("E") || line.starts_with("N") && !line.contains("E")
    {
    fastq_2_sequence.push(line);
   } else if line.starts_with("+") || line.starts_with("-") {
       fastq_2_strand.push(line);
   } else if line.contains("E"){
       fastq_2_quality.push(line);
   }

  }

   let mut fastq_1_pre:Vec<FileFastqPre> = Vec::new();
   let mut fastq_2_pre:Vec<FileFastqPre> = Vec::new();

   for i in 0..fastq_1_header.len(){
       fastq_1_pre.push(FileFastqPre{
           header: fastq_1_header[i].clone(),
           sequence: fastq_1_sequence[i].clone(),
           strand: fastq_1_strand[i].clone(),
           quality:fastq_1_quality[i].clone(),
       })
   }

   for i in 0..fastq_2_header.len()
   {
      fastq_2_pre.push(FileFastqPre{
          header: fastq_2_header[i].clone(),
          sequence:fastq_2_sequence[i].clone(),
          strand:fastq_2_strand[i].clone(),
          quality:fastq_2_quality[i].clone(),
      })
   }

   #[derive(Debug,Clone)]
   struct FileFastqPost {
       header: String,
       sequence:String,
       strand: String,
       quality:String,
   }

   let fastq1_clean:Vec<FileFastqPost> = Vec::new();
    for i in fastq_1.iter_mut(){
        let mut tupleseq = i.sequence.clone().chars().collect::<Vec<_>>();
        let mut tuplequality = i.quality.clone().chars().collect::<Vec<_>>();
            let mut tupleclean = Vec::new();
            let mut tuplequalityclean = Vec::new();
            for j in 0..tupleseq.len() {
                  if tuplequality[j] == '@' {
                  continue
                } else {
                    tupleclean.push(tupleseq[j].to_string());
                    tuplequalityclean.push(tuplequality[j].to_string());
                }
        }
        println!("{:?}, {:?}", tupleclean.join(""), tuplequalityclean.join(""));
    }

   }
