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
	
	data.frame(
		row.names=attr(seq, 'name'),
		mut=mut,
		pot.mut=length(pot.mut),
		ctrl=ctrl,
		pot.ctrl=length(pot.ctrl),
		p.value=fisher.test(m,alternative='greater')$p.value
	)
}

op <- OptionParser(
	usage=
"
RHypermut: hypermutation detection.
Detects APOBEC3-induced hypernutation in HIV sequences. Works the same as the LANL Hypermut progran (https://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html).

usage: %prog --fasta=FASTA [--cleanfasta=CLEANFASTA] [--hyperfasta=HYPERFASTA] [--log=LOG] [--makeconsensus] [--verbose]"
)
op <- add_option(
	op,
	"--fasta",
	help="Input fasta file. If --makeconsensus is set the first sequence is used as a reference otherwise a concensus sequence is created to use as a reference.",
	type='character'
)
op <- add_option(
	op,
	"--cleanfasta",
	help="Output fasta file with suspected hypermutated sequences removed.",
	type='character',
	default=NA
)
op <- add_option(
	op,
	"--hyperfasta",
	help="Output fasta file containing the suspected hypermutated sequences. (Optional).",
	type='character',
	default=NA
)
op <- add_option(
	op,
	"--log",
	help="CSV file that contains the number of hypermutations, potential hypermutations, control mutations, potential control mutations and p-value for each sequence. (Optional).",
	type='character',
	default=NA
)
op <- add_option(
	op,
	"--makeconsensus",
	help="Set to use make a consensus sequence to use as a reference. (Optional).",
	type='logical',
	action='store_true',
	default=F
)
op <- add_option(
	op,
	"--verbose",
	help="Set to print the log file to screen. (Optional).",
	type='logical',
	action='store_true',
	default=F
)
args <- parse_args(op)

fasta.file <- args$fasta
clean.fasta.file <- args$cleanfasta
hyper.fasta.file <- args$hyperfasta
log.file <- args$log
make.consensus <- args$makeconsensus
verbose <- args$verbose

f <- lapply(read.fasta(fasta.file), tolower)

if (make.consensus) {
	ref <- consensus(do.call(rbind, f))
} else {
	ref <- f[[1]]
	f <- f[-1]
}

data <- do.call(rbind, lapply(f, check.hyper, ref=ref))

p.value <- data$p.value

if (verbose)
	print(data)

if (!is.na(log.file))
	write.csv(data, log.file)

if (!is.na(clean.fasta.file)) {
	f.clean <- f[p.value >= 0.05] 
	write.fasta(f.clean, names(f.clean), clean.fasta.file)
}

if (!is.na(hyper.fasta.file)) {
	f.hyper <- f[p.value < 0.05]
	write.fasta(f.hyper, names(f.hyper), hyper.fasta.file)
}
