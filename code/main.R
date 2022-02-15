library('seqinr') # Biological Sequences Retrieval and Analysis
library('digest') # Compact Hash Digests of R Objects
source('code/functions.R')

r6 <- read.fasta(file = 'data/R6.fa')

genome = r6$AE007317.1
genome = paste(toupper(genome), collapse='')

### Reading k-mers --------------
k = 14

length(genome)

genome <- 'TTGAAAGAAAAACAATTTTG'
kmer(genome, k=14)

### Simple Jaccard distances --------------
A = c('TTGAAAGAAAAACA','TGAAAGAAAAACAA','GAAAGAAAAACAAT',
      'AAAGAAAAACAATT','AAGAAAAACAATTT','AGAAAAACAATTTT','GAAAAACAATTTTG')
B = c('TTGAAAGAAAAACA','TGAAAGAAAAACAA','GAAAGAAAAACAAT',
      'AAAGAAAAACAATT','GAAAAACAATTTTG','CTCGATCCATGTAT','TCGATCCATGTATG')

length(intersect(A,B))/length(union(A,B))
length(intersect(A,B))
length(union(A,B))

### MinHash Jaccard distances --------------
# objective: take an ubbiased sample of the k-mers and calculate Jaccard distance from this smaller subset
hash = unlist(lapply(A, digest, algo = "murmur32", serialize = F, seed = 0))
sort(hash)

s = 1000 # sketch size