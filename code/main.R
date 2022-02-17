# load necessary packages
library('seqinr')  # Biological Sequences Retrieval and Analysis
library('seqR')    # Fast and Comprehensive K-Mer Counting Package
library('digest')  # Compact Hash Digests of R Objects
library('ape')     # Analyses of Phylogenetics and Evolution
library('dplyr')   # data manipulation
library('tidyr')   # tidy messy data
library('ggplot2') # data viz
#library('kmer')   # Fast K-Mer Counting
source('code/functions.R') # load own functions

k = 14 # size of kmer

# Read and process fasta into kmers ---------------------
# fasta2kmer is defined in 'code/functions.R'
R6_kmers = fasta2kmer(file = 'data/R6.fa', k)
TIGR4_kmers = fasta2kmer(file = 'data/TIGR4.fa', k)

# name 'draft' sequences as D82, D84
D82_kmers = fasta2kmer(file = 'data/14412_3#82.contigs_velvet.fa', k)
D84_kmers = fasta2kmer(file = 'data/14412_3#84.contigs_velvet.fa', k)

### create hash and sketches from kmers ---------------------
# create vector of input names
inputs = c('R6', 'TIGR4', 'D82' ,'D84')
inputs = factor(inputs, levels=inputs)
ninputs = length(inputs) # number of inputs

# create sorted hashes for all inputs
for (i in inputs){
  print(paste('input:', i)) 
  
  # create hash and sort
  assign(paste(i,'hash',sep='_'),
         sort(unlist(lapply(get(paste(i,'kmers',sep='_')), 
                       digest, algo = "murmur32", serialize = F, seed = 0)), method='quick'))
}

ss = c(1E3, 1E4, 1E5, 1E6) # varying sketch sizes

# create sketches with varying sketch sizes for all inputs
for (s in ss){
  print(paste('sketch size:', format(s, scientific=T)))
  
  for (i in inputs){
    print(paste('input:', i)) 
    
    # extract sketch from hash
    sketch_name = paste(i,'sketch',format(s, scientific=T),sep='_')
    assign(sketch_name,
           get(paste(i,'hash',sep='_'))[1:s])
    
    # save sketch in .txt file
    # cat(get(sketch_name),
    #    file=paste0('output/',i,'_sketch_',format(s, scientific=T),'.txt'), sep="\n")
  }
}

### compute MinHash Jaccard distances --------------
# for varying sketch sizes for all inputs
res_sketch = list() # to store results
r = 1 # row counter
for (s in ss){
  #! calculation is doubled here, not ideal
  for (i in 1:ninputs){
    for (j in 1:ninputs){
      # calculate Jaccard distance and store results in a row
      res_sketch[[r]] = c(format(s, scientific=T), 
                       as.character(inputs[i]), 
                       as.character(inputs[j]),
                       J_dist(get(paste(inputs[i],'sketch',format(s, scientific=T),sep='_')),
                              get(paste(inputs[j],'sketch',format(s, scientific=T),sep='_'))))
      r = r+1
    }
  }
}

res_sketch = as.data.frame(do.call(rbind, res_sketch)) # convert list to data frame

# from complete hash
res_hash = list() # to store results
r = 1 # row counter
for (i in 1:ninputs){
  for (j in 1:ninputs){
    # calculate Jaccard distance and store results in a row
    res_hash[[r]] = c("H", 
                     as.character(inputs[i]), 
                     as.character(inputs[j]),
                     J_dist(get(paste(inputs[i],'hash',sep='_')),
                            get(paste(inputs[j],'hash',sep='_'))))
    r = r+1
  }
}

# convert list to data frame
res_hash = as.data.frame(do.call(rbind, res_hash))

# combine results into single data frame
res = rbind(res_hash, res_sketch)
names(res) = c('s', 'x', 'y', 'J') # apply column names
res = res %>% mutate(J = as.numeric(J)) # convert distance to numeric

### Plot of distances ----------------
# differences by sample
p1 <- ggplot(mapping = aes(x=factor(y, levels=inputs), y=J, col=s)) +
  # plot J distance for varying s
  geom_point(data = filter(res, s != 'H', J > 0)) +
  # plot J distance computed from whole hash
  geom_point(data = filter(res, s == 'H', J >0), pch=4, col=1, size=3) +
  facet_wrap(~factor(x, levels=inputs), scales='free_y') +
  labs(title = 'Comparison of Jaccard distance between samples for varying s',
       subtitle = 'Full distances are marked by crosses') +
  ylab('Jaccard distance') +
  xlab('Sample') +
  theme_bw() + 
  theme(legend.position = 'bottom')

# differences by s
p2 <- ggplot(mapping = aes(x=s, y=J)) +
  # plot J distance for varying s
  geom_point(data = filter(res, s != 'H', J > 0)) +
  # plot J distance computed from whole hash
  geom_point(data = filter(res, s == 'H', J > 0), pch=4, col=1, size=3) +
  facet_wrap(~factor(x, levels=inputs) + factor(y, levels=inputs), scales='free_y') +
  labs(title = 'Comparison of Jaccard distance between samples for varying s',
       subtitle = 'Full distances are marked by crosses') +
  ylab('Jaccard distance') +
  xlab('Sample') +
  theme_bw() + 
  theme(legend.position = 'bottom')

ggsave(filename = 'output/figures/J_dist_sample.png', p1, width = 8, height = 6)
ggsave(filename = 'output/figures/J_dist_s.png', p2, width = 10, height = 6)

### Neighbour Joining tree ---------------------
# convert distances in data frame to matrix
J_mats = list()
for (ss in unique(res$s)){
  J_mats[[ss]] = res %>%
    filter(s == ss) %>%
    pivot_wider(-s, names_from = y, values_from = J)
}
J_mats = lapply(J_mats, df.matrix) 

# plot neighbour joining trees
png(file='output/figures/nj_hash.png',
    width=600, height=400)
plot(ape::nj(J_mats$H))
dev.off()

png(file='output/figures/nj_sketch.png',
    width=600, height=400)
plot(ape::nj(J_mats$`1e+03`))
dev.off()

