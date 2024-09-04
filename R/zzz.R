.onLoad <- function(...) {
  cmdstan_version <- cmdstanr::cmdstan_version(error_on_NA = FALSE)
  if (is.null(cmdstan_version)) {
    stop(paste0("No CmdStan installation found.",
                "Run cmdstanr::install_cmdstan() to install."),
         call. = FALSE)
  }
}
## From https://github.com/rok-cesnovar/misc/blob/master/democmdstanr/R/zzz.R
