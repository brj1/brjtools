library(seqinr)
library(optparse)

#implementation of Hypermut (https://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html)
get.pot.mut.nuc <- function(ref, seq) {
	l <- length(seq)
	
	which(
		seq[2:(l-1)] %in% c('a', 'g') &
			seq[3:l] != 'c' &
			ref[1:(l-2)] == 'g'
	)
}

get.pot.ctrl.nuc <- function(ref, seq) {
	l <- length(seq)
	
	which(
		((seq[2:(l-1)] %in% c('a', 'g') & seq[3:l] == 'c') |
			seq[2:(l-1)] %in% c('c', 't')) &
			ref[1:(l-2)] == 'g'
	)
}

get.pot.mut.iupac <- function(ref, seq) {
  l <- length(seq)
  
  which(
    (seq[2:(l-1)] %in% c('a', 'g') | seq[2:(l-1)] == 'r') &
      seq[3:l] != 'c' &
      ref[1:(l-2)] == 'g'
  )
}

get.pot.ctrl.iupac <- function(ref, seq) {
  l <- length(seq)
  
  which(
    (((seq[2:(l-1)] %in% c('a', 'g') | seq[2:(l-1)] == 'r') & seq[3:l] == 'c') |
       (seq[2:(l-1)] %in% c('c', 't') | seq[2:(l-1)] == 'y')) &
      ref[1:(l-2)] == 'g'
  )
}

get.pot.mut <- get.pot.mut.iupac
get.pot.ctrl <- get.pot.ctrl.iupac

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
		rate.ratio=(mut / length(pot.mut)) / (ctrl / length(pot.ctrl)),
		p.value=fisher.test(m,alternative='greater')$p.value
	)
}

try.open <- function(file.name, FOO, file.type, ...) {
	tryCatch(
		FOO(file.name, ...),
		error=function() {
			error(
				sprintf(
					"%s file '%s' does not exist or cannot be opened.",
					file.type,
					file.name
				)
			)
		}
	)
}

try.write <- function(x, file.name, FOO, file.type, ...) {
	tryCatch(
		FOO(x, file.name, ...),
		error=function() {
			error(
				sprintf(
					"Cannot create or write to %s file '%s'.",
					file.type,
					file.name
				)
			)
		}
	)
}

op <- OptionParser(
	usage=
"\r            
RHypermut: hypermutation detection.
Detects APOBEC3-induced hypernutation in HIV sequences. Works the same as the LANL Hypermut progran (https://www.hiv.lanl.gov/content/sequence/HYPERMUT/hypermut.html).

Usage: %prog [--cleanfasta=CLEANFASTA] [--hyperfasta=HYPERFASTA] [--log=LOG] [--ref=REF] [--threshold=THRESHOLD] [--makeconsensus] [--quiet] FASTA

\tFASTA
\t\tInput fasta file."
)
op <- add_option(
	op,
	"--cleanfasta",
	help="Output fasta file with suspected hypermutated sequences removed. Does not include the reference sequence.",
	type='character',
	default=NA
)
op <- add_option(
	op,
	"--hyperfasta",
	help="Output fasta file containing the suspected hypermutated sequences.",
	type='character',
	default=NA
)
op <- add_option(
	op,
	"--log",
	help="CSV file that contains the number of hypermutations, potential hypermutations, control mutations, potential control mutations and p-value for each sequence.",
	type='character',
	default=NA
)
op <- add_option(
	op,
	"--ref",
	help="Either taxon name of the reference sequence in FASTA or a fasta file containing a single reference sequence. If --makeconcensus is specified this is not used. If not specified then the first sequence of FASTA is used as the reference. Ignored if --makeconsensus is set.",
	type='character',
	default='1'
)
op <- add_option(
	op,
	"--threshold",
	help="The p-value threshold to assign sequences as hypermutated. Default: 0.05.",
	type='numeric',
	default='0.05'
)
op <- add_option(
	op,
	"--makeconsensus",
	help="Set to make a consensus sequence to use as a reference.",
	type='logical',
	action='store_true',
	default=F
)
op <- add_option(
	op,
	"--quiet",
	help="Set to suppress output.",
	type='logical',
	action='store_true',
	default=F
)
all.args <- parse_args(op, positional_arguments=1)
args <- all.args$options

fasta.file <- all.args$args
clean.fasta.file <- args$cleanfasta
hyper.fasta.file <- args$hyperfasta
log.file <- args$log
ref.name <- args$ref
threshold <- args$threshold
make.consensus <- args$makeconsensus
quiet <- args$quiet

f <- lapply(try.open(fasta.file, read.fasta, 'input fasta'), tolower)

if (make.consensus) {
	ref <- consensus(do.call(rbind, f))
} else {
	suppressWarnings(ref.num <- as.integer(ref.name))
	if (!is.na(ref.num) && ref.num <= length(names(f)) && ref.num > 0) {
		names(f)[ref.num] <- ref.name
	}
	
	if (ref.name %in% names(f)) {
		ref <- f[[ref.name]]
		f <- f[names(f) != ref.name]
	} else {
		ref <- tolower(
			try.open(ref.name, read.fasta, 'reference fasta')[[1]]
		)
	}
}

if (!all(sapply(f, length) == length(ref))) {
	error("Sequence lengths differ. Make sure that the input fasta file is aligned and that the reference sequence is the same length.")
}

data <- do.call(rbind, lapply(f, check.hyper, ref=ref))
p.value <- data$p.value

if (!quiet)
	print(data)

if (!is.na(log.file))
	try.write(data, log.file, write.csv, 'log')

if (!is.na(clean.fasta.file)) {
	f.clean <- f[p.value >= threshold] 
	try.write(
		f.clean,
		names=names(f.clean),
		clean.fasta.file,
		write.fasta,
		'clean output fasta'
	)
}

if (!is.na(hyper.fasta.file)) {
	f.hyper <- f[p.value < threshold]
	try.write(
		f.hyper,
		names=names(f.hyper),
		hyper.fasta.file,
		write.fasta,
		'hypermutated output fasta'
	)
}
