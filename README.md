# PhyloAcc-C

Genomes contain conserved non-coding elements (CNEs) that are important to biological function. But which elements relate to which traits, and how? PhyloAcc-C answers this question by estimating how molecular substitution rates and the rate of change of a continuous trait covary.

Input to PhyloAcc-C is a multiple sequence alignment of a CNE, continuous trait data observed in extant species, and a background phylogeny and substitution process. Gibbs sampling is used to assign latent conservation states (background, conserved, accelerated) to lineages. These hidden states are associated with molecular and phenotypic change by rate multipliers. Covariation in rates of molecular and trait evolution is assessed using the posterior ratio of trait variance multipliers.

For examples of how to use the PhyloAcc-C R package, see the [tutorial](tutorial/tutorial.pdf).

For further information on the PhyloAcc-C model, see the [model description](model_description.pdf).

Once you install the PhyloAcc-C package in R, help is available in the usual way, e.g. `?pac3` brings up a jump page to documentation on each function.

## Installation instructions

PhyloAcc-C has been implemented using a mixture of R and C++. This is necessary to get good performance but has the downside that you will need a development environment installed on your computer. (For the time being, at least. In future we will provide binary packages.)

First make sure you are set up for using R and C++.

Windows: install [Rtools](https://cran.r-project.org/bin/windows/Rtools/index.html).

Mac: install Xcode by typing the command `xcode-select --install` on Terminal.

Ubuntu: use command `sudo apt-get install r-base-dev`.

Now download the [source R package](pac3_0.1.tar.gz). Open R studio and click on the "Packages" tab. Click "Install". Click "Browse..." and select the source package you just downloaded and then chose "Install". R should build and install the PhyloAcc-C package.

## Further information

More information on the PhyloAcc family of methods is available at the [PhyloAcc organization](https://phyloacc.github.io).
