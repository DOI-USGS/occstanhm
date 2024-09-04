## FROM artifactory.wma.chs.usgs.gov/docker-official-mirror/rocker/verse
FROM rocker/verse

# COPY DOIRootCA2.crt /usr/local/share/ca-certificates
# RUN update-ca-certificates
# ENV CURL_CA_BUNDLE /etc/ssl/certs/ca-certificates.crt
# RUN SSL_CERT_FILE=~/trusted-certs.pem

## Adapted from https://hub.docker.com/r/asachet/rocker-stan/dockerfil2
## To strip down the rocker version
## and also write the install_rstan.R file inside of docker instead of
## copying it locally
RUN apt-get update \
	&& apt-get install -y --no-install-recommends \
	clang \
	gdebi

## From https://github.com/stan-dev/rstan/wiki/Installing-RStan-on-Linux
RUN echo \
	'dotR <- file.path(Sys.getenv("HOME"), ".R")\n\
	if (!file.exists(dotR)) dir.create(dotR)\n\
	M <- file.path(dotR, "Makevars")\n\
	if (!file.exists(M)) file.create(M)\n\
	cat("\nCXX14FLAGS=-O3 -march=native -mtune=native -fPIC",\n\
	"CXX14=g++",\n\
	file = M, sep = "\n", append = TRUE)'\
	| cat >> install_rstan.R

RUN ["r", "install_rstan.R"]


# Install rstan
RUN install2.r --error --deps TRUE \
    rstan \
    && rm -rf /tmp/downloaded_packages/ /tmp/*.rds

# Installing the rest
RUN install2.r -s --error \
	bayesplot \
	brms \
	coda \
	loo \
	projpred \
	rstanarm \
	rstantools \ 
	shinystan \
	tidybayes \
	ggthemes \
	GGally \
	ggpubr \
	data.table \
	sf \
	betareg \
	patchwork \
	ggthemes \
	kableExtra \
	flextable \
	lintr \
  bibtex \
  knitcitations

RUN R -e "update.packages(repos = 'http://cran.us.r-project.org')"

RUN apt install -y gcc-10 g++-10
RUN  update-alternatives --install /usr/bin/gcc gcc /usr/bin/gcc-10 100 --slave /usr/bin/g++ g++ /usr/bin/g++-10 --slave /usr/bin/gcov gcov /usr/bin/gcov-10
RUN update-alternatives --config gcc

RUN R -e 'remove.packages(c("rstan"))'
RUN R -e 'remove.packages(c("StanHeaders"))'
RUN R -e 'remotes::install_github("stan-dev/rstan/StanHeaders@experimental", upgrade = "always", force = TRUE)'
RUN R -e 'remotes::install_github("stan-dev/rstan/rstan/rstan@experimental", upgrade = "always", force = TRUE)'

RUN wget --progress=dot:mega https://github.com/stan-dev/cmdstan/releases/download/v2.34.1/cmdstan-2.34.1.tar.gz
RUN tar -zxpf cmdstan-2.34.1.tar.gz
RUN ln -s cmdstan-2.34.1 cmdstan
RUN cd cmdstan; make stan-update; make build -j4


RUN R -e 'remotes::install_github("stan-dev/cmdstanr")'
RUN R -e 'cmdstanr::install_cmdstan()'

RUN curl -LO https://quarto.org/download/latest/quarto-linux-amd64.deb
RUN gdebi --non-interactive quarto-linux-amd64.deb

RUN CURL_CA_BUNDLEA=~/trusted-certs.pem
