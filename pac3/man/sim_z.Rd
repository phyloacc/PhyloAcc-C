\name{sim_z}
\alias{sim_z}
\title{
Simulate latent conservation states under the PhyloAcc-C model
}
\description{
Generates latent conservation states. Requires specification of a tree,
transition probabilities, and a categorical prior on the conservation state of the
branch leading to the root node.

Can be used to explore and visualize your prior on latent conservation states. Its
return value is a necessary input for simulating nucleotide alignments or traits.
}
\usage{
sim_z(tree, Z_R_prior, P_nc, P_ca, P_ac, plot_it=F)
}
\arguments{
\item{tree}{A single tree as read by the \code{\link[ape]{read.tree}} function.}
\item{Z_R_prior}{Specifies the probability the branch leading to the root node is in state 1, 2, or 3. Use a length 3 numeric vector.}
\item{P_nc}{The probability a branch becomes conserved given a neutral/background predecessor.}
\item{P_ca}{The probability a branch becomes accelerated given a conserved predecessor.}
\item{P_ac}{The probability a branch reverts to a conserved state given an accelerated predecessor. Set to 0 for a Dollo process.}
\item{plot_it}{When TRUE, display tree with edges coloured by generated conservation state (black=1/bg/neutral, blue=2/conserved, red=3/accelerated).}
}
\value{
A vector specifying conservation states for branches that can be passed to \code{\link{sim_X}} or \code{\link{sim_y}}.
}
\examples{
x<-1+1 # Todo
}
