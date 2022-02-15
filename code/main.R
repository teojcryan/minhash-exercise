library('seqinr') # Biological Sequences Retrieval and Analysis
library('seqR') # Fast and Comprehensive K-Mer Counting Package
library('digest') # Compact Hash Digests of R Objects
#library('kmer')
source('code/functions.R')

### Read data -------------
R6 <- read.fasta(file = 'data/R6.fa')
R6_genome = paste(toupper(R6[[1]]), collapse='') # extract genome as string

TIGR4 <- read.fasta(file = 'data/TIGR4.fa')
TIGR4_genome = paste(toupper(TIGR4[[1]]), collapse='') # extract genome as string

### Reading k-mers --------------
k = 14

R6_kmer <- kmer2(R6_genome, k)
TIGR4_kmer <- kmer2(TIGR4_genome, k)

# Attempts using own function/ other packages
# r6_kmer <- kmer(genome_string, k) # own function takes too long
# r6_kmer <- seqinr::count(genome, k) # Error: memory exhausted (limit reached?)
# r6_kmer <- kmer::kcount(genome, k) # Error: cannot allocate vector of size 14.0 Gb

### Simple Jaccard distances --------------
A = c('TTGAAAGAAAAACA','TGAAAGAAAAACAA','GAAAGAAAAACAAT',
      'AAAGAAAAACAATT','AAGAAAAACAATTT','AGAAAAACAATTTT','GAAAAACAATTTTG')
B = c('TTGAAAGAAAAACA','TGAAAGAAAAACAA','GAAAGAAAAACAAT',
      'AAAGAAAAACAATT','GAAAAACAATTTTG','CTCGATCCATGTAT','TCGATCCATGTATG')

J_dist(A,B)

### MinHash Jaccard distances --------------
# objective: take an unbiased sample of the k-mers and calculate Jaccard distance from this smaller subset
s = 10000 # sketch size

R6_hash = unlist(lapply(R6_kmer, digest, algo = "murmur32", serialize = F, seed = 0))
R6_sketch = sort(R6_hash, method="quick")[1:s]
cat(R6_sketch,file="output/R6_sketch.txt",sep="\n")

TIGR4_hash = unlist(lapply(TIGR4_kmer, digest, algo = "murmur32", serialize = F, seed = 0))
TIGR4_sketch = sort(TIGR4_hash, method="quick")[1:s]
cat(TIGR4_sketch,file="output/TIGR4_sketch",sep="\n")

J_dist(R6_hash, TIGR4_hash)
