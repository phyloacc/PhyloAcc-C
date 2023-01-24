\name{sim_X}
\alias{sim_X}
\title{
Simulate a nucleotide alignment under the PhyloAcc-C model
}
\description{
Generates an alignment of length S sites. Requires specification of a
nucleotide substitution process, a tree, latent conservation states, and
nucleotide rate multipliers.
}
\usage{
sim_X(PI, Q, tree, neu_branches, con_branches, acc_branches, con_m=0.1, acc_m=0.7, S=200)
}
\arguments{
\item{PI}{The stationary distribution of the nucleotide substitution process, used to pick ancestral nucleotides at the root node. Provide a 4 element numeric vector.}
\item{Q}{A 4x4 numeric rate matrix specifying the nucleotide substitution process. Index order is a c g t.}
\item{tree}{A single tree as read by the \code{\link[ape]{read.tree}} function.}
\item{neu_branches}{Character vector of node names identifying background/neutral branches.}
\item{con_branches}{Character vector of node names identifying conserved branches.}
\item{acc_branches}{Character vector of node names identifying accelerated branches.}
\item{con_m}{The conserved substitution rate multiplier i.e. \eqn{r_2}.}
\item{acc_m}{The accelerated substitution rate multiplier i.e. \eqn{r_3}.}
\item{S}{Number of sites to simulate i.e. the length of the alignment.}
}
\value{
A nucleotide alignment in character matrix form ('a','c','g','t') that can be passed to \code{\link{run_model}}.
}
\seealso{
Use \code{\link{save_fasta}} if you wish to save the result of \code{sim_X} for viewing in other software.
}
\examples{
x<-1+1 # Todo
}
