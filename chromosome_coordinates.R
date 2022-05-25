if (!require("BiocManager", quietly = TRUE))
  install.packages("BiocManager")

BiocManager::install("ensembldb")
BiocManager::install('AnnotationHub')

library(ensembldb)
library(AnnotationHub)

ah <- AnnotationHub()
ensdb = query(ah, c("EnsDb", "v104", "Homo sapiens"))

df <- read.csv('first_table.txt', sep='\t')

second_table = data.frame('gene_name'=character(), 'domain'=character(), 'enst_id'=character(), 
                          'ensp_id'=character(), 'start'=character(), 'end'=character(), 
                          'part_type'=character(), 'source'=character(), 'chromosome'=integer(),'start_chr'=integer(), 'end_chr'=integer())

data = df
data$start_chr = 25
data$end_chr = 25
db = ensdb[[1]]
for(i in 1:nrow(data)){
  x <- data[i,]
  chr = x$chromosome
  edbx <- filter(db, filter = ~ seq_name == as.character(chr))
  glob_start = c()
  glob_end = c()
  prot_coordinates = length(strsplit(x$start,',')[[1]])
  for (i in c(1:prot_coordinates)){
    start = strtoi(strsplit(x$start,',')[[1]][i])
    end = strtoi(strsplit(x$end,',')[[1]][i])
    name = x$ensp_id
    gene = IRanges(start = start, end = end, names = name)
    output = proteinToGenome(gene, edbx)
    chr_start  = output[[1]]@ranges@start
    chr_width = output[[1]]@ranges@width
    chr_end = c()
    chromosome_coordinates = length(chr_start)
    for (j in c(1:chromosome_coordinates)){
      n = chr_start[j] + chr_width[j] - 1
      chr_end = c(chr_end,n)
    }
    glob_start = c(glob_start, chr_start)
    glob_end = c(glob_end, chr_end)
  }
  x$start_chr = paste(glob_start, collapse = ', ')
  x$end_chr = paste(glob_end, collapse = ', ')
  print(x)
  second_table <- rbind(second_table,x)
}

write.table(second_table, file="second_table.txt", sep = '\t', row.names = FALSE)