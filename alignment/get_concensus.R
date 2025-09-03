#!/usr/bin/Rscript
library(seqinr)
library(optparse)
library(dplyr)

get.mode <- function(x, dup=NA) {
  tab <- table(x)
  m <- max(tab)
  if (!is.na(dup) & sum(tab == m) > 1)
    dup
  else
    names(tab)[tab == m]
}

args <- OptionParser() %>%
  add_option("--name", default = "CONCENSUS") %>%
  parse_args2()

fasta.file <- args$args[1]
fasta.concensus.file <- args$args[2]
name <- args$options$name

f <- read.fasta(fasta.file)
df <- do.call(rbind, f)
cons <- list(apply(df, 2, get.mode, '-'))

write.fasta(cons, name, fasta.concensus.file)
