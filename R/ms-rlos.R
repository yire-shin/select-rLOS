# Top-level interface for the msrlos package.

#' Model selection for r-largest order statistics
#'
#' Runs the model-selection procedures for r-largest order-statistic data.
#'
#' @param xdat Numeric matrix. Each row is a block (for example, a year), and
#'   columns contain the largest observations in decreasing order.
#' @param sig.ed Significance level used by the entropy-difference procedures.
#' @param dmin Minimum spacing used when constructing fixed h candidates.
#' @param numh Number of candidate h values considered when h is selected
#'   automatically.
#' @return An object of class `msrlos` with components `surv`, `rmed`,
#'   `alt.surv`, and `alt.rmed`.
#' @author Shin Yire, Jihong Park, and Jeong-Soo Park
#'   (correspondence: jspark@jnu.ac.kr).
#' @export
ms.rlos <- function(xdat, sig.ed = 0.1, dmin = 0.03, numh = 8) {
  if (is.data.frame(xdat)) xdat <- as.matrix(xdat)
  if (!is.matrix(xdat) || !is.numeric(xdat)) {
    stop("`xdat` must be a numeric matrix or data frame.", call. = FALSE)
  }
  if (nrow(xdat) < 2L || ncol(xdat) < 2L) {
    stop("`xdat` must have at least two rows and two order-statistic columns.", call. = FALSE)
  }
  if (!is.numeric(sig.ed) || length(sig.ed) != 1L || is.na(sig.ed) || sig.ed <= 0 || sig.ed >= 1) {
    stop("`sig.ed` must be a single number between 0 and 1.", call. = FALSE)
  }
  if (!is.numeric(dmin) || length(dmin) != 1L || is.na(dmin) || dmin < 0) {
    stop("`dmin` must be a single non-negative number.", call. = FALSE)
  }
  if (!is.numeric(numh) || length(numh) != 1L || is.na(numh) || numh < 2) {
    stop("`numh` must be at least 2.", call. = FALSE)
  }
  numh <- as.integer(numh)

  R <- ncol(xdat)

  batch <- redtest.all(
    xdat, h.fix = NULL, sig.ed = sig.ed,
    dmin = dmin, numh = numh
  )

  surv <- sel.surv(
    xdat, sig.ed = sig.ed, h.fix = batch$h.fix,
    dmin = dmin, survtest = batch
  )

  rmed <- sel.rmed(
    xdat, sig.ed = sig.ed, h.fix = batch$h.fix,
    dmin = dmin, medtest = batch
  )

  alt.surv <- Alteria.new(
    xdat, sig.ed = sig.ed, start = "surv",
    mid.best = surv$mstar, rhat.best = surv$rstar,
    h.fix = batch$h.fix, altest = batch, altering = TRUE
  )

  otr <- alt.surv$rstar
  bid <- alt.surv$mstar
  rcol <- if (otr == 1) R else otr - 1
  alt.surv$theta <- batch$redtest[[bid]][rcol, c(4:6, 10)]

  alt.rmed <- Alteria.new(
    xdat, sig.ed = sig.ed, start = "rmed",
    mid.best = rmed$mstar, rhat.best = rmed$rstar,
    h.fix = batch$h.fix, altest = batch, altering = TRUE
  )

  otr <- alt.rmed$rstar
  bid <- alt.rmed$mstar
  rcol <- if (otr == 1) R else otr - 1
  alt.rmed$theta <- batch$redtest[[bid]][rcol, c(4:6, 10)]

  z <- list(
    surv = surv,
    rmed = rmed,
    alt.surv = alt.surv,
    alt.rmed = alt.rmed
  )

  class(z) <- c("msrlos", "list")
  z
}

#' Print msrlos results
#'
#' @param x An object returned by [ms.rlos()].
#' @param ... Unused.
#' @return `x`, invisibly.
#' @export
print.msrlos <- function(x, ...) {
  pick <- function(obj) {
    if (is.null(obj)) {
      return(c(model = NA_character_, r = NA_character_))
    }
    model <- if (!is.null(obj$model.names) && !is.null(obj$mstar) &&
                 length(obj$model.names) >= obj$mstar) {
      obj$model.names[obj$mstar]
    } else if (!is.null(obj$best.model)) {
      sub("\\s+[0-9]+$", "", as.character(obj$best.model)[1])
    } else if (!is.null(obj$model.names)) {
      obj$model.names[1]
    } else {
      NA_character_
    }
    r <- if (!is.null(obj$rstar)) as.character(obj$rstar) else NA_character_
    c(model = as.character(model)[1], r = r[1])
  }

  out <- rbind(
    Survival = pick(x$surv),
    `r-Median` = pick(x$rmed),
    `Alter Surv` = pick(x$alt.surv),
    `Alter r-Med` = pick(x$alt.rmed)
  )

  cat("msrlos model-selection results\n\n")
  print(out, quote = FALSE)
  invisible(x)
}
