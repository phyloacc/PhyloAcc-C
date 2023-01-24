\name{pac3-package}
\alias{pac3-package}
\alias{pac3}
\docType{package}
\title{
PhyloAcc-C method for relating molecular and phenotypic evolution
}
\description{

\packageDescription{pac3}

Further information about the PhyloAcc family of models is available at \url{https://phyloacc.github.io}.

}
\details{

Given a tree and substitution process, immediately apply PhyloAcc-C to alignment
and continuous trait data using \code{\link{run_model}}, then examine the trace using
\code{\link{plot_mixing}}. Plot results using \code{\link{plot_v}} and \code{\link{plot_z}}.

Preferably, first simulate representative data (using \code{\link{sim_X}},
\code{\link{sim_y}}, and \code{\link{sim_z}}) to see if the model works well on your tree.

}
\author{
\packageAuthor{pac3}

Maintainer: \packageMaintainer{pac3}
}
\references{

PhyloAcc-C: Gemmell, P., Sackton, T. B., Edwards, S. V., and Liu, J. S. (YYYY) A phylogenetic method linking nucleotide substitution rates to continuous phenotypic innovation, \emph{in preparation}.

Original Phylo-Acc paper: Hu, Z., Sackton, T. B., Edwards, S. V., and Liu, J. S. (2019) Bayesian detection of convergent rate changes of conserved noncoding elements on phylogenetic trees, \emph{Molecular biology and evolution}, 36(5), 1086--1100.

}
\keyword{package}
\seealso{
Use \code{\link[ape:ape-package]{ape}} to manipulate trees (esp. \code{\link[ape]{read.tree}}) and alignments (esp. \code{\link[ape]{read.dna}}) before applying PhyloAcc-C.
}
