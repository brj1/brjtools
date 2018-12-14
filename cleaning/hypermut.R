library(seqinr)
library(optparse)

#implementation of Hypermut (https://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html)
get.pot.mut <- function(ref, seq) {
	l <- length(seq)
	
	which(
		seq[2:(l-1)] %in% c('a', 'g') &
			seq[3:l] != 'c' &
			ref[1:(l-2)] == 'g'
	)
}

get.pot.ctrl <- function(ref, seq) {
	l <- length(seq)
	
	which(
		((seq[2:(l-1)] %in% c('a', 'g') & seq[3:l] == 'c') |
			seq[2:(l-1)] %in% c('c', 't')) &
			ref[1:(l-2)] == 'g'
	)
}

check.hyper <- function(ref, seq) {
	pot.mut <- get.pot.mut(ref, seq)
	pot.ctrl <- get.pot.ctrl(ref, seq)
	mut <- sum(seq[pot.mut] == 'a')
	ctrl <- sum(seq[pot.ctrl] == 'a')
	
	m <- t(
		matrix(
			c(mut, length(pot.mut) - mut, ctrl, length(pot.ctrl) - ctrl),
			nrow=2
		)
	)
	
	if (verbose)
		cat(mut, length(pot.mut), ctrl, length(pot.ctrl), '\n')
	
	fisher.test(m, alternative='greater')$p.value
}

op <- OptionParser()
op <- add_option(op, "--fasta", type='character')
op <- add_option(op, "--cleanfasta", type='character')
op <- add_option(op, "--hyperfasta", type='character', default=NA)
op <- add_option(op, "--makeconsensus", type='logical', action='store_true', default=F)
op <- add_option(op, "--verbose", type='logical', action='store_true', default=F)
args <- parse_args(op)

fasta.file <- args$fasta
clean.fasta.file <- args$cleanfasta
hyper.fasta.file <- args$hyperfasta
make.consensus <- args$makeconsensus
verbose <- args$verbose

f <- lapply(read.fasta(fasta.file), tolower)

if (make.consensus) {
	ref <- consensus(do.call(rbind, f))
} else {
	ref <- f[[1]]
	f <- f[-1]
}

p.value <- lapply(f, check.hyper, ref=ref)

if (verbose)
	print(data.frame(row.names=names(p.value), p.value=unlist(p.value)))

f.clean <- f[p.value >= 0.05] 
write.fasta(f.clean, names(f.clean), clean.fasta.file)

if (!is.na(hyper.fasta.file)) {
	f.hyper <- f[p.value < 0.05]
	write.fasta(f.hyper, names(f.hyper), hyper.fasta.file)
}