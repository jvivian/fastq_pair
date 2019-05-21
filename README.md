# Fastq-Pairing Algorithms in Rust

Attempting to benchmark different fastq-pairing methods to identify trade-offs between memory and runtime.


| Method         | Description                                                                               | Source                       |
|----------------|-------------------------------------------------------------------------------------------|------------------------------|
| Store-Read     | Iter R1 and store. Iter R2 and write out.                                                 | https://tinyurl.com/ya2tcg8k |
| Iter-Both      | Iter through both R1/R2 simultaneously.  Write/pop hashmap as pairs are found.            | https://tinyurl.com/ya7l3amo |
| Seek-Read      | Hash headers to byte position.  Iter R2 and seek to R1 to write out pairs.                | Joel                         |
| Seek-Iter-Both | Iter through both R1/R2 storing byte position. Seek/Write/Pop hashmap as pairs are found. | Joel/John                    |


## Implemented
- [x] Store-Read
- [x] Iter-Both
- [x] Seek-Read
- [ ] Seek-Iter-Both

## Benchmarked
Here are benchmarks on a pair of ~6GB fastq files where one had 10% of reads shuffled relative to the second file

| Method     | Time (m) | Memory (G) |
|------------|----------|------------|
| EdwardsLab | 9.9      | 16         |
| Store      | 2.1      | 7.5        |
| Seek       | 4.25     | 1.9        |
| Iter       | 1.45     | 0.0        |


## Additional features
- [ ] BAM input
- [x] GZIP input / output
- [ ] Assert paired end
- [x] Singletons
- [ ] Include non-unique header descriptors
- [ ] Derived/custom output names
