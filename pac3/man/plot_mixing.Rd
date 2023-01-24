\name{plot_mixing}
\alias{plot_mixing}
\title{
Plot trace of nucleotide and variance multipliers from MCMC output
}
\description{
Plot trace of log(\eqn{beta_3}) - log(\eqn{beta_2}) as well as rate multipliers 
(\eqn{r_2}, \eqn{r_3}), entries of Phi matrix 
(c, a, b), and variance multipliers (\eqn{sigma^2}, \eqn{log(beta_2)}, \eqn{log(beta_3)}).
}
\usage{
plot_mixing(mcmc_out, ith=1)
}
\arguments{
\item{mcmc_out}{The list produced after calling \code{\link{run_model}}.}
\item{ith}{Plot every ith value, useful when \code{mcmc_out} contains many samples.}
}
\examples{
\dontrun{
# run the PhyloAcc-C model for 5000 iterations
mcmc_out = pac3::run_model(verbose=T, tree=tree, alignment=alignment, trait_tbl=trait_tbl,
                           transition_matrix=Q, stationary_distro=PI, n_iter=5000,
                           inner_iter=500, Phi_prior=c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0))

# plot a trace
plot_mixing(mcmc_out)

# plot a trace, thinning to every 5th draw
plot_mixing(mcmc_out, ith=5)
}
}