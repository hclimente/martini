#' Get shake settings.
#' 
#' @description Creates a list composed by all \code{shake} settings, with the values provided by the user, or the default ones if none is provided.
#' @param associationScore Association score to measure association between genotype and phenotype. Possible values: skat (default), chi2, trend.
#' @param modelScore Model selection criterion Possible values: consistency, bic (default), aic, aicc, mbic.
#' @param encoding SNP encoding. Possible values: additive (default), resessive, dominant, codominant.
#' @param debug Display additional information. Possible values: TRUE, FALSE (default).
#' @return A list of \code{shake} settings.
get_evo_settings <- function(...){

  settings <- list(...)

  # unsigned int
  if (! "associationScore" %in% names(settings)) {
    settings[["associationScore"]] = 0
  } else if (settings[["associationScore"]] == "skat") {
    settings[["associationScore"]] = 0
  } else if (settings[["associationScore"]] == "chi2") {
    settings[["associationScore"]] = 1
  } else if (settings[["associationScore"]] == "trend") {
    settings[["associationScore"]] = 2
  } else {
    stop(paste("Error: invalid associationScore", settings[["associationScore"]]))
  }

  # unsigned int
  if (! "modelScore" %in% names(settings)) {
    settings[["modelScore"]] = 1
  } else if (settings[["modelScore"]] == "consistency") {
    settings[["modelScore"]] = 0
  } else if (settings[["modelScore"]] == "bic") {
    settings[["modelScore"]] = 1
  } else if (settings[["modelScore"]] == "aic") {
    settings[["modelScore"]] = 2
  } else if (settings[["modelScore"]] == "aicc") {
    settings[["modelScore"]] = 3
  } else if (settings[["modelScore"]] == "mbic") {
    settings[["modelScore"]] = 4
  } else {
    stop(paste("Error: invalid modelScore", settings[["modelScore"]]))
  }
  
  # unsigned int
  if (! "encoding" %in% names(settings)) {
    settings[["encoding"]] = 0;
  } else if (settings[["encoding"]] == "additive") {
    settings[["encoding"]] = 0
  } else if (settings[["encoding"]] == "recessive") {
    settings[["encoding"]] = 1
  } else if (settings[["encoding"]] == "dominant") {
    settings[["encoding"]] = 2
  } else if (settings[["encoding"]] == "codominant") {
    settings[["encoding"]] = 3
  } else {
    stop(paste("Error: invalid encoding", settings[["encoding"]]))
  }
  
  # bool
  if (! "debug" %in% names(settings)) {
    settings[["debug"]] = FALSE;
  } else if (! is.logical(settings[["debug"]])) {
    stop("Error: debug must be logical.")
  }

  return(settings);
}
