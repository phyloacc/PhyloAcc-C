\name{plot_v}
\alias{plot_v}
\title{
Plot posterior distribution of log(beta_3) - log(beta_2)
}
\description{
Display the posterior distribution of the relative rate of phenotypic change for
accelerated versus conserved branches as a histogram.
}
\usage{
plot_v(mcmc_out, drop_n=1000, thin_n=1)
}
\arguments{
\item{mcmc_out}{The list produced after calling \code{\link{run_model}}.}
\item{drop_n}{Drop the first n elements as burn-in.}
\item{thin_n}{Use only every nth element, useful when \code{mcmc_out} contains many samples.}
}
\examples{
\dontrun{
# run the PhyloAcc-C model for 5000 iterations
mcmc_out = pac3::run_model(verbose=T, tree=tree, alignment=alignment, trait_tbl=trait_tbl,
                           transition_matrix=Q, stationary_distro=PI, n_iter=5000,
                           inner_iter=500, Phi_prior=c(1.0, 1.0, 1.0, 1.0, 1.0, 1.0))

# plot histogram of log b3 - log b2, discarding 1000 elements as burn-in.
plot_v(mcmc_out, drop_n=1000)
}
}