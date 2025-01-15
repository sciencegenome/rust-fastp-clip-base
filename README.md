# rust-fastp-clip-base
 - rust-fastp-clip-base
 - specifically dropping quality at specific bases.
 - It generates 4 files ones with the dropped quality at that base and the ones with the all quality dropped below that base. 
 - please see the last commit message and if it says compiled binary then it is completed or else still in development version.

 ```
 cargo build 

 ```
 ```
 ╭─gauravsablok@fedora ~/Downloads/rust-fastp-clip-base-main  
 ╰─➤  ./rust-fastp-clip-base -h
Usage: rust-fastp-clip-base <READS_1_ARG> <READS_2_ARG> <QUALITY_SCORE>

Arguments:
  <READS_1_ARG>    please provide the reads R1 file path
  <READS_2_ARG>    please provide the reads R2 file path
  <QUALITY_SCORE>  please provide the quality value to be used as a threshold

Options:
  -h, --help     Print help
  -V, --version  Print version
 - to run the compiled binary. 
 
 ```
  ./target/debug/rust-fastp-clip-base ./sample-files/test1.fastq ./sample-files/test2.fastq 30
 
 ```
 
 Gaurav Sablok
