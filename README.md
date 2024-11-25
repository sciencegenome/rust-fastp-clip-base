# rust-fastp-clip-base
 - rust-fastp-clip-base
 - specifically dropping quality at specific bases.
 - It generates 4 files ones with the dropped quality at that base and the ones with the all quality dropped below that base. 

 ```
 cargo build 

 ```
  ./target/debug/rust-fastp-clip-base ./sample-files/test1.fastq ./sample-files/test2.fastq 30
 
 ```
 Gaurav Sablok
