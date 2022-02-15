kmer <- function(genome, k, toList = F){
  hash <- new.env(hash = TRUE)
  
  for (i in 1:(nchar(genome)-k+1)){
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

J_dist <- function(A, B){
  # Compute the Jaccard distance of two input samples A, B
  1 - length(intersect(A, B)) / length(union(A, B))
}