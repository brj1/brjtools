library(seqinr)
library(tidyverse)

ggmsa <- function(fasta.file, consensus=NA) {
f <- read.fasta(fasta.file)

f.matrix <- do.call(rbind, f)

if (!is.na(consensus)) {
  if (consensus %in% row.names(f.matrix)) {
    cons <- f.matrix[consensus, ]
  } else {
    cons <- consensus
  }
  
  f.matrix <- apply(
    f.matrix,
    1,
    function(x) {
      x[tolower(x) == tolower(cons)] <- NA
      x
    }
  ) %>%
  	t
}

f.matrix <- as.data.frame(f.matrix)
f.matrix$name <- row.names(f.matrix)

f.fort <- gather(f.matrix, "position", "cod", -name)

ggplot(f.fort, aes(x=position, y=name, fill=cod)) +
	geom_tile(show.legend=F) +
	scale_fill_discrete(na.value=NA) +
	theme_bw() +
	theme(panel.grid.major.x=element_blank(), panel.grid.minor.x=element_blank()) +
	scale_x_continuous(name="Codon Position", expand=c(0, 0)) +
	scale_y_discrete(name="Sequence ID")
}
