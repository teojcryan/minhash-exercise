library('seqinr') # Biological Sequences Retrieval and Analysis
library('digest') # Compact Hash Digests of R Objects
library('kmer')
source('code/functions.R')

r6 <- read.fasta(file = 'data/R6.fa')

genome = r6

### Reading k-mers --------------
k = 14
genome_string = paste(toupper(genome[[1]]), collapse='')
r6_kmer <- kmer(genome_string, k) # takes too long

# Try methods from existing packages
# r6_kmer <- seqinr.count(genome, k) # Error: memory exhausted (limit reached?)
# r6_kmer <- kmer.kcount(genome, k) # Error: cannot allocate vector of size 14.0 Gb

### Simple Jaccard distances --------------===
A = c('TTGAAAGAAAAACA','TGAAAGAAAAACAA','GAAAGAAAAACAAT',
      'AAAGAAAAACAATT','AAGAAAAACAATTT','AGAAAAACAATTTT','GAAAAACAATTTTG')
B = c('TTGAAAGAAAAACA','TGAAAGAAAAACAA','GAAAGAAAAACAAT',
      'AAAGAAAAACAATT','GAAAAACAATTTTG','CTCGATCCATGTAT','TCGATCCATGTATG')

J_dist(A,B)

### MinHash Jaccard distances --------------
# objective: take an unbiased sample of the k-mers and calculate Jaccard distance from this smaller subset
s = 5 # sketch size

hash_A = unlist(lapply(A, digest, algo = "murmur32", serialize = F, seed = 0))
sketch_A = sort(hash_A)[1:s]
hash_B = unlist(lapply(B, digest, algo = "murmur32", serialize = F, seed = 0))
sketch_B = sort(hash_B)[1:s]

sketch_A
sketch_B
J_dist(sketch_A, sketch_B)
