# rust-fastp-clip-base
 - rust-fastp-clip-base
 - specifically dropping quality at specific bases.
 - It generates 4 files ones with the dropped quality at that base and the ones with the all quality dropped below that base. 
 
 ```
 cargo build 

 ```
 - to run the compiled binary. 
 
 ```
  ./target/debug/rust-fastp-clip-base ./sample-files/test1.fastq ./sample-files/test2.fastq 30
 
 ```
 - general note: Incase of Golang and RUST, please see the last commit message and if it says compiled binary then it is completed or else still in development version.  
 
 Gaurav Sablok
