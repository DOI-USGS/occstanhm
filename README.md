# occstanhm: Hierarchical occupancy models with correlated error structure

#### Authors:          Richard A. Erickson, Charles J. Labuzzetta
#### Point of contact: Richard A. Erickson (rerickson@usgs.gov)
#### Repository Type:  _R_ packing calling _Stan_ models
#### Year of Origin:   2024 (original publication)
#### Year of Version:  2024
#### Version:          2.0.0 
#### Digital Object Identifier (DOI): XXXXXX
#### USGS Information Product Data System (IPDS) no.: IP-XXXXX (internal agency tracking)

***

_Suggested Citation:_

Erickson, RA, and Labuzzetta, CJ.
2024.
`occstanhm`: Hierarchical occupancy models with correlated error structure.
U.S. Geological Survey software release.
Reston, Va.
https://doi.org/10.5066/XXXXXXX.

_Authors' [ORCID](https://orcid.org) nos.:_

- Richard A. Erickson, [0000-0003-4649-482X](https://orcid.org/0000-0003-4649-482X)
- Charles J. Labuzzetta, [0000-0002-6027-0120](https://orcid.org/0000-0002-6027-0120)

***
***

This repository contains a R package with 2-level and 3-level occupancy models.
The package also contains a set of tutorials designed to help users understand how to use and code occupancy models in Stan.
The models include a statistical hierarchy that allows for
multi-species modeling and estimating correlations among species.
The models may also be used as a more general "random-effect" type
occupancy model.

The models are written in Stan and called through R.
The `cmdstanr` package is used, rather than the more common `rstan`
package, because features used in Stan were not supported by `rstan`
at the time of development.
Additionally, `cmdstanr` allows for quicker performance than `rstan`.

# Installation

This code requires the `cmdstanr` package to run in R.
Please look up the [official documentation](https://mc-stan.org/cmdstanr/) for [install directions](https://mc-stan.org/cmdstanr/articles/cmdstanr.html).
We also include a `Dockerfile` for people who prefer to use Docker. 

Once you have `cmdstanr` installed and setup, this package may be installed using this code:

```{r}
if (!require("remotes")) install.packages("remotes")
remotes::install_gitlab('umesc/quant-ecology/occstanhm@main', host='code.usgs.gov')
```

You may wish to lockdown a specific version by changing `main` to the version (e.g., `v2.0`).
As of May 2024, the `build_vignettes` option appears to not be working.
Vigenttes may be built by cloning the repository using git and then installing locally.

# Where to get started

After installing the program, please consult the vignettes to learn more.
The `Introduction_overview` provides a starting place and describes a suggested learning path through the tutorials.

# Repository Files

This repository contains the code for an R package using RStan.
This repository contains the standard R repository files
(see the [official R Documentation Writing R Extensions](https://cran.r-project.org/doc/manuals/r-release/R-exts.html)
accessed May 2024 or the online book, [_R Pakcages (2e)_](https://r-pkgs.org/) for an descriptions of these
files).
In addition to the R Package source files, this repository contains
the following files:

This repository file contains the following files and folder:

- `README.md` is this file.
- `LICENSE.md` is the Official USGS License. 
- `code.json` is the code metadata.
- `CONTRIBUTING.md` describes how to contribute to this project.
- `DISCLAIMER.md` is the standard USGS disclaimer.
- `.gitignore` is a file telling git which files to not track.
- `docker_files` contains the `Dockerfile` to use this code.


# User skill level

This package expects a user to understand Bayesian Statistics and occupanyc models.
Users seeking to adapt code would also benfit from understanding programming in Stan.

# Acknowledgments

This research was funded by the USGS Biological Threats and 
Invasive Species Research Program and the US Fish and
Wildlife Service.
Any use of trade, firm, or product names is for descriptive purposes only
and does not imply endorsement by the U.S. Government.


[r_ext]: https://cran.r-project.org/doc/manuals/r-release/R-exts.html
[stan_user_1.13]: https://mc-stan.org/docs/2_22/stan-users-guide/multivariate-hierarchical-priors-section.html
[DF_header]: https://code.usgs.gov/-/ide/project/rerickson/climatchr/edit/master/-/README.md#dockerfile-contents
[DF_build]: https://code.usgs.gov/-/ide/project/rerickson/climatchr/edit/master/-/README.md#build-dockerfile
[DF_run]: https://code.usgs.gov/-/ide/project/rerickson/climatchr/edit/master/-/README.md#docker-run
