### Utility functions

fasta2kmer <- function(filename, k){
  # loads a fasta file and reads all relevant kmers
  fasta <- read.fasta(file = filename)
  kmers = unique(unlist(lapply(fasta, function(x){
    genome = paste(toupper(x), collapse='')
    kmer2(genome, k)
  })))
  return(kmers)
}

kmer2 <- function(genome_string, k){
  # reads kmers from a genome string
  kmers <- seqR::count_kmers(genome_string, k,
                               with_kmer_counts=TRUE, # to track counts of kmers (not essential)
                               with_kmer_names=TRUE)
  kmers <- colnames(kmers) # extract column names
  kmers <- gsub('[[:punct:][:digit:]]', '', kmers) # make kmers human readable
  return(kmers)
}

kmer <- function(genome_string, k, toList = F){
  # Original function to convert genome string to kmers, but too slow, use kmer2 instead
  hash <- new.env(hash = TRUE)
  
  for (i in 1:(nchar(genome_string)-k+1)){
    # extract relevant kmer
    kmer <- substr(genome, i, i+k-1)
    
    # if kmer does not exist in hash, add key with new value = 1 (first count)
    if (!(kmer %in% ls(hash))) {
      hash[[kmer]] <- 1
    } else {
    # if kmer-key exists in hash, add 1 to value
    hash[[kmer]] <- hash[[kmer]] + 1
    }
    
    if (i %% 1E5 == 0) {paste('iteration', i)}
  }
  if (toList == T){as.list(hash)} else {hash}
}

# Attempts using own function/ other packages
# kmer(genome_string, k) # own function takes too long
# seqinr::count(genome_vec, k) # Error: memory exhausted (limit reached?)
# kmer::kcount(genome_vec, k) # Error: cannot allocate vector of size 14.0 Gb

J_dist <- function(A, B){
  # Compute the Jaccard distance of two input samples A, B
  1 - length(intersect(A, B)) / length(union(A, B))
}

df.matrix <- function(df){
  # convert a data frame into a matrix
  m <- as.matrix(df[-1])
  rownames(m) <- t(df[,1]) # convert first column into matrix rownames
  m
}
