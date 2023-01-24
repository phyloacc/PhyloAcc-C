\name{save_fasta}
\alias{save_fasta}
\title{
Save a simulated alignment to a FASTA file
}
\description{
Convenience function to save results of sim_X in FASTA file format.
}
\usage{
save_fasta(X_char, file)
}
\arguments{
\item{X_char}{A character matrix output by sim_X.}
\item{file}{The name of a file to write to.}
}
\examples{
\dontrun{
# save the alignment X
save_fasta(X, "X.fna")
}
}
