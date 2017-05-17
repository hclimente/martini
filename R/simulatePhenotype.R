simulatePhenotype <- function(X, snps, h2, model = "additive", effectSize = rnorm(sum(snps)), qualitative = FALSE, ncases, ncontrols){
  # check correspondence with gcta implementation

  X <- X[, snps]
  
  # get effect sizes u and weights w
  u <- effectSize
  p <- (2 * colSums(X == 2) + colSums(X == 1)) / (2 * nrow(X))
  x <- 2 * (X == 2) + (X == 1)
  w = (x - 2 * p) / sqrt(2 * p * (1 - p))
  
  if (model == "additive")
    geno <- colSums(t(w) * u)
  else
    stop(paste0("Genetic model ", model, " not recognised."))
  
  residual.var <- var(geno) * (1 / h2 - 1)
  residual <- rnorm(length(geno), sd = sqrt(residual.var))
  
  trait <- geno + residual
  
  if (qualitative){
    if (! exists("ncases") )
      stop("Specify ncases if qualitative = TRUE.")
    else if (! exists("ncontrols") )
      stop("Specify ncontrols if qualitative = TRUE.")
    else if ( length(trait) < (ncases + ncontrols) )
      stop("Cases and controls requested exceed number of samples provided.")

    trait.sorted <- sort(trait, index.return = TRUE)
    cases <- head(trait.sorted$ix, n = ncases)
    controls <- tail(trait.sorted$ix, n = ncontrols)
    
    Y <- numeric(length(trait))
    Y <- NA
    Y[cases] <- 1
    Y[controls] <- 0
  } else {
    Y <- trait
  }
  
  return(Y)
}