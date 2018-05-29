library(seqinr)

args <- commandArgs(trailingOnly = T)

fasta.file <- args[1]
fasta.dedup.file <- args[2]

f <- read.fasta(fasta.file)

cat("Fasta read:",  fasta.file, "\n")

same <- lapply(f, function(x) which(unlist(lapply(f, function(y) all(x == y)))))
same <- same[unlist(lapply(same, length)) > 1]
same <- unique(lapply(same, function(x) {names(x) <- NULL; x}))

cat("\nFound", length(same), "sets of duplicate sequences:\n")
sup <- lapply(same, function(x) cat(names(f)[x], '\n'))

dont.keep <- unlist(lapply(same, function(x) x[-1]))

f.dedup <- f[-dont.keep]

cat("\n", length(dont.keep), " sequences removed:\n", sep="")
cat(names(f)[dont.keep], "\n\n")

write.fasta(f.dedup, names(f.dedup), fasta.dedup.file)

cat("Unique seqences written to:", fasta.dedup.file, "\n")