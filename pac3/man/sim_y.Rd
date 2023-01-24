\name{sim_y}
\alias{sim_y}
\title{
Simulate a continuous trait under the PhyloAcc-C model
}
\description{
Generates trait values. Requires specification of a tree, latent conservation states, 
variance multipliers, and the ancestral value of the trait at the root node.
}
\usage{
sim_y(tree, neu_branches, con_branches, acc_branches, v=1, b2=0.2, b3=2, yR=0)
}
\arguments{
\item{tree}{A single tree as read by the \code{\link[ape]{read.tree}} function.}
\item{neu_branches}{Character vector of node names identifying background/neutral branches.}
\item{con_branches}{Character vector of node names identifying conserved branches.}
\item{acc_branches}{Character vector of node names identifying accelerated branches.}
\item{v}{Baseline \eqn{\sigma^2} for trait movement.}
\item{b2}{The multiplier \eqn{\beta_2} applied to \eqn{\sigma^2} on conserved branches.}
\item{b3}{The multiplier \eqn{\beta_3} applied to \eqn{\sigma^2} on accelerated branches.}
\item{yR}{Ancestral value of the trait at the root node.}
}
\value{
A data frame containing the value of the simulated trait at each node of the tree; can be passed to \code{\link{run_model}}.
}
\examples{
x<-1+1 # Todo
}
