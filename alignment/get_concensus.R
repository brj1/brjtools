#!/usr/bin/Rscript
library(seqinr)

get.mode <- function(x, dup=NA) {
  tab <- table(x)
  m <- max(tab)
  if (!is.na(dup) & sum(tab == m) > 1)
    dup
  else
    names(tab)[tab == m]
}

args <- commandArgs(trailingOnly = T)

fasta.file <- args[1]
fasta.concensus.file <- args[2]

f <- read.fasta(fasta.file)
df <- do.call(rbind, f)
cons <- list(apply(df, 2, get.mode, '-'))

write.fasta(cons, "CONCENSUS", fasta.concensus.file)
