hurdle <- function(formula, data, subset, na.action,
                   dist = c("poisson", "negbin", "geometric"),
                   zero.dist = c("binomial", "poisson", "negbin", "geometric"),
                   link = c("logit", "probit", "cloglog", "cauchit", "log"),
		   control = hurdle.control(...),
		   model = TRUE, y = TRUE, x = FALSE, ...)
{
  ## set up likelihood components
  zeroPoisson <- function(parms) {
    ## mean
    lambda <- as.vector(exp(Z %*% parms))
    ## log-likelihood
    loglik0 <- -lambda ## = dpois(0, lambda, log = TRUE)
    ## collect and return
    loglik <- sum(loglik0[Y0]) + sum(log(1 - exp(loglik0[Y1])))
    loglik
  }

  countPoisson <- function(parms) {
    ## mean
    lambda <- as.vector(exp(X %*% parms))[Y1]
    ## log-likelihood
    loglik0 <- -lambda ## = dpois(0, lambda, log = TRUE)
    loglik1 <- dpois(Y[Y1], lambda, log = TRUE)
    ## collect and return
    loglik <- sum(loglik1) - sum(log(1 - exp(loglik0)))
    loglik
  }

  zeroNegBin <- function(parms) {
    ## parameters
    lambda <- as.vector(exp(Z %*% parms[1:kz]))
    theta <- exp(parms[kz+1])
    ## log-likelihood
    loglik0 <- suppressWarnings(dnbinom(0, size = theta, mu = lambda, log = TRUE))
    ## collect and return
    loglik <- sum(loglik0[Y0]) + sum(log(1 - exp(loglik0[Y1])))
    loglik
  }

  countNegBin <- function(parms) {
    ## parameters
    lambda <- as.vector(exp(X %*% parms[1:kx]))[Y1]
    theta <- exp(parms[kx+1])
    ## log-likelihood
    loglik0 <- suppressWarnings(dnbinom(0, size = theta, mu = lambda, log = TRUE))
    loglik1 <- suppressWarnings(dnbinom(Y[Y1], size = theta, mu = lambda, log = TRUE))
    ## collect and return
    loglik <- sum(loglik1) - sum(log(1 - exp(loglik0)))
    loglik
  }

  zeroGeom <- function(parms) zeroNegBin(c(parms, 0))
  countGeom <- function(parms) countNegBin(c(parms, 0))

  zeroBinom <- function(parms) {
    ## mean
    lambda <- as.vector(linkinv(Z %*% parms))
    ## log-likelihood
    loglik <- sum(log(1 - lambda[Y0])) + sum(log(lambda[Y1]))
    loglik  
  }

  dist <- match.arg(dist)
  zero.dist <- match.arg(zero.dist)
  countDist <- switch(dist,
                     "poisson" = countPoisson,
		     "geometric" = countGeom,
		     "negbin" = countNegBin)
  zeroDist <- switch(zero.dist,
                     "poisson" = zeroPoisson,
		     "geometric" = zeroGeom,
		     "negbin" = zeroNegBin,
		     "binomial" = zeroBinom)
  loglikfun <- function(parms) countDist(parms[1:(kx + (dist == "negbin"))]) +
    zeroDist(parms[(kx + (dist == "negbin") + 1):(kx + kz + (dist == "negbin") + (zero.dist == "negbin"))])

  ## binary link processing
  link <- match.arg(link)
  linkinv <- make.link(link)$linkinv

  if(control$trace) cat("Hurdle Count Model\n",
    paste("count model:", dist, "with log link\n"),
    paste("zero hurdle model:", zero.dist, "with", ifelse(zero.dist == "binomial", link, "log"), "link\n"),
    sep = "")


  ## call and formula
  cl <- match.call()
  if(missing(data)) data <- environment(formula)
  mf <- match.call(expand.dots = FALSE)
  m <- match(c("formula", "data", "subset", "na.action"), names(mf), 0)
  mf <- mf[c(1, m)]
  mf$drop.unused.levels <- TRUE
  
  ## extended formula processing
  if(length(formula[[3]]) > 1 && identical(formula[[3]][[1]], as.name("|")))
  {
    ff <- formula
    formula[[3]][1] <- call("+")
    mf$formula <- formula
    ffc <- . ~ .
    ffz <- ~ .
    ffc[[2]] <- ff[[2]]
    ffc[[3]] <- ff[[3]][[2]]
    ffz[[3]] <- ff[[3]][[3]]
    ffz[[2]] <- NULL
  } else {
    ffz <- ffc <- ff <- formula
    ffz[[2]] <- NULL
  }
  if(any(sapply(unlist(as.list(ffz[[2]])), function(x) identical(x, as.name("."))))) {
    ffz <- eval(parse(text = sprintf( paste("%s -", deparse(ffc[[2]])), deparse(ffz) )))
  }

  ## call model.frame()
  mf[[1]] <- as.name("model.frame")
  mf <- eval(mf, parent.frame())
  
  ## extract terms, model matrices, response
  mt <- terms(formula, data = data)
  mtX <- terms(ffc, data = data)
  X <- model.matrix(mtX, mf)
  mtZ <- terms(ffz, data = data)
  mtZ <- terms(update(mtZ, ~ .), data = data)
  Z <- model.matrix(mtZ, mf)
  Y <- model.response(mf, "numeric")


  ## sanity checks
  if(length(Y) < 1) stop("empty model")
  if(all(Y > 0)) stop("invalid dependent variable, minimum count is not zero")  
  if(!isTRUE(all.equal(as.vector(Y), as.integer(round(Y + 0.001)))))
    stop("invalid dependent variable, non-integer values")
  Y <- as.integer(round(Y + 0.001))
  if(any(Y < 0)) stop("invalid dependent variable, negative counts")
  
  if(control$trace) {
    cat("dependent variable:\n")
    tab <- table(factor(Y, levels = 0:max(Y)), exclude = NULL)
    names(dimnames(tab)) <- NULL
    print(tab)
  }

  ## convenience variables
  n <- length(Y)
  kx <- NCOL(X)
  kz <- NCOL(Z)
  Y0 <- Y < 1
  Y1 <- Y > 0
    
  ## starting values
  start <- control$start
  if(!is.null(start)) {
    valid <- TRUE
    if(!("count" %in% names(start))) {
      valid <- FALSE
      warning("invalid starting values, count model coefficients not specified")
      start$count <- rep(0, kx)
    }
    if(!("zero" %in% names(start))) {
      valid <- FALSE
      warning("invalid starting values, zero-inflation model coefficients not specified")
      start$zero <- rep(0, kz)
    }
    if(length(start$count) != kx) {
      valid <- FALSE
      warning("invalid starting values, wrong number of count model coefficients")
    }
    if(length(start$zero) != kz) {
      valid <- FALSE
      warning("invalid starting values, wrong number of zero-inflation model coefficients")
    }
    if(dist == "negbin" | zero.dist == "negbin") {
      if(!("theta" %in% names(start))) start$theta <- c(1, 1)
      start <- list(count = start$count, zero = start$zero, theta = rep(start$theta, length.out = 2))
      if(is.null(names(start$theta))) names(start$theta) <- c("count", "zero")      
      if(dist != "negbin") start$theta <- start$theta["zero"]
      if(zero.dist != "negbin") start$theta <- start$theta["count"]
    } else {
      start <- list(count = start$count, zero = start$zero)
    }
    if(!valid) start <- NULL
  }
  
  if(is.null(start)) {
    if(control$trace) cat("generating starting values...")
    model_count <- glm.fit(X, Y, family = poisson())
    model_zero <- switch(zero.dist,
      "poisson" = glm.fit(Z, Y, family = poisson()),
      "negbin" = glm.fit(Z, Y, family = poisson()),
      "geometric" = suppressWarnings(glm.fit(Z, factor(Y > 0), family = binomial())),
      "binomial" = suppressWarnings(glm.fit(Z, factor(Y > 0), family = binomial(link = link))))
    start <- list(count = model_count$coefficients, zero = model_zero$coefficients)
    start$theta <- c(count = if(dist == "negbin") 1 else NULL,
                     zero = if(zero.dist == "negbin") 1 else NULL)
    if(control$trace) cat("done\n")
  }

  ## model fitting
  ## control parameters
  method <- control$method
  hessian <- control$hessian
  separate <- control$separate
  ocontrol <- control
  control$method <- control$hessian <- control$separate <- control$start <- NULL

  ## ML estimation
  ## separate estimation of censored and truncated component...
  if(separate) {
    if(control$trace) cat("calling optim() for count component estimation:\n")
    fit_count <- optim(fn = countDist,
      par = c(start$count, if(dist == "negbin") log(start$theta["count"]) else NULL),
      method = method, hessian = hessian, control = control)

    if(control$trace) cat("calling optim() for zero hurdle component estimation:\n")
    fit_zero <- optim(fn = zeroDist,
      par = c(start$zero,  if(zero.dist == "negbin") log(start$theta["zero"]) else NULL),
      method = method, hessian = hessian, control = control)
    if(control$trace) cat("done\n")
    fit <- list(count = fit_count, zero = fit_zero)

    ## coefficients
    coefc <- fit_count$par[1:kx]
    coefz <- fit_zero$par[1:kz]
    theta <- c(count = if(dist == "negbin") as.vector(exp(fit_count$par[kx+1])) else NULL,
               zero = if(zero.dist == "negbin") as.vector(exp(fit_zero$par[kz+1])) else NULL)
    ## covariances
    vc_count <- -solve(as.matrix(fit_count$hessian))
    vc_zero <- -solve(as.matrix(fit_zero$hessian))
    SE.logtheta <- list()
    if(dist == "negbin") {
      SE.logtheta$count <- as.vector(sqrt(diag(vc_count)[kx+1]))
      vc_count <- vc_count[-(kx+1), -(kx+1), drop = FALSE]
    }
    if(zero.dist == "negbin") {
      SE.logtheta$zero <- as.vector(sqrt(diag(vc_zero)[kz+1]))
      vc_zero <- vc_zero[-(kz+1), -(kz+1), drop = FALSE]
    }
    vc <- rbind(cbind(vc_count, matrix(0, kx, kz)), cbind(matrix(0, kz, kx), vc_zero))
    SE.logtheta <- unlist(SE.logtheta)
  } else {
  ## ...or joint.
    if(control$trace) cat("calling optim() for joint count and zero hurlde estimation:\n")
    fit <- optim(fn = loglikfun,
      par = c(start$count, if(dist == "negbin") log(start$theta["count"]) else NULL,
              start$zero,  if(zero.dist == "negbin") log(start$theta["zero"]) else NULL),
      method = method, hessian = hessian, control = control)
    if(fit$convergence > 0) warning("optimization failed to converge")
    if(control$trace) cat("done\n")

    ## coefficients
    coefc <- fit$par[1:kx]
    coefz <- fit$par[(kx + (dist == "negbin") + 1):(kx + kz + (dist == "negbin"))]
    ## covariances
    vc <- -solve(as.matrix(fit$hessian))
    np <- c(if(dist == "negbin") kx+1 else NULL,
            if(zero.dist == "negbin") kx+kz+1+(dist == "negbin") else NULL)
    if(length(np) > 0) {
      theta <- as.vector(exp(fit$par[np]))
      SE.logtheta <- as.vector(sqrt(diag(vc)[np]))
      names(theta) <- names(SE.logtheta) <- c(if(dist == "negbin") "count" else NULL,
                                              if(zero.dist == "negbin") "zero" else NULL)
      vc <- vc[-np, -np, drop = FALSE]
    } else {
      theta <- NULL
      SE.logtheta <- NULL
    }

  }
  names(coefc) <- names(start$count) <- colnames(X)
  names(coefz) <- names(start$zero) <- colnames(Z)
  colnames(vc) <- rownames(vc) <- c(paste("count", colnames(X), sep = "_"),
                                    paste("zero",  colnames(Z), sep = "_"))

  ## fitted and residuals
  phi <- if(zero.dist == "binomial") linkinv(Z %*% coefz)[,1] else exp(Z %*% coefz)[,1]
  p0_zero <- switch(zero.dist,
                    "binomial" = 1-phi,
		    "poisson" = exp(-phi),
		    "negbin" = dnbinom(0, size = theta["zero"], mu = phi),
		    "geometric" = dnbinom(0, size = 1, mu = phi))

  mu <- exp(X %*% coefc)[,1]
  p0_count <- switch(dist,
		    "poisson" = exp(-mu),
		    "negbin" = dnbinom(0, size = theta["count"], mu = mu),
		    "geometric" = dnbinom(0, size = 1, mu = mu))

  Yhat <- (1 - p0_zero) * mu / (1 - p0_count)
  res <- Y - Yhat


  rval <- list(coefficients = list(count = coefc, zero = coefz),
    residuals = res,
    fitted.values = Yhat,
    optim = fit,
    method = method,
    control = control,
    start = start,
    n = n,
    df.null = n - 2,
    df.residual = n - (kx + kz + 2 * (dist == "negbin")),
    terms = list(count = mtX, zero = mtZ, full = mt),
    theta = theta,
    SE.logtheta = SE.logtheta,
    loglik = if(separate) fit_count$value + fit_zero$value else fit$value,
    vcov = vc,
    dist = list(count = dist, zero = zero.dist),
    link = if(zero.dist == "binomial") link else NULL,
    linkinv = if(zero.dist == "binomial") linkinv else NULL,
    separate = separate,
    converged = if(separate) fit_count$convergence < 1 & fit_zero$convergence < 1 else fit$convergence < 1,
    call = cl,
    formula = ff,
    levels = .getXlevels(mt, mf),
    contrasts = list(count = attr(X, "contrasts"), zero = attr(Z, "contrasts"))
  )
  if(model) rval$model <- mf
  if(y) rval$y <- Y
  if(x) rval$x <- list(count = X, zero = Z)
      
  class(rval) <- "hurdle"
  return(rval)
}

hurdle.control <- function(method = "BFGS", maxit = 10000, trace = FALSE, separate = TRUE, start = NULL, ...) {
  rval <- list(method = method, maxit = maxit, trace = trace, separate = separate, start = start)
  rval <- c(rval, list(...))
  if(!is.null(rval$fnscale)) warning("fnscale must not be modified")
  rval$fnscale <- -1
  if(!is.null(rval$hessian)) warning("hessian must not be modified")
  rval$hessian <- TRUE
  if(is.null(rval$reltol)) rval$reltol <- .Machine$double.eps^(1/1.6)
  rval
}

coef.hurdle <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  rval <- object$coefficients
  rval <- switch(model,
                 "full" = structure(c(rval$count, rval$zero),
		   .Names = c(paste("count", names(rval$count), sep = "_"),
                   paste("zero", names(rval$zero), sep = "_"))),
		 "count" = rval$count,
		 "zero" = rval$zero)
  rval
}

vcov.hurdle <- function(object, model = c("full", "count", "zero"), ...) {
  model <- match.arg(model)
  rval <- object$vcov
  if(model == "full") return(rval)

  cf <- object$coefficients[[model]]
  wi <- seq(along = object$coefficients$count)
  rval <- if(model == "count") rval[wi, wi] else rval[-wi, -wi]
  colnames(rval) <- rownames(rval) <- names(cf)
  return(rval)
}

logLik.hurdle <- function(object, ...) {
  structure(object$loglik, df = object$n - object$df.residual, class = "logLik")
}

print.hurdle <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "\n", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    cat(paste("Count model coefficients (truncated ", x$dist$count, " with log link):\n", sep = ""))
    print.default(format(x$coefficients$count, digits = digits), print.gap = 2, quote = FALSE)
    if(x$dist$count == "negbin") cat(paste("Theta =", round(x$theta["count"], digits), "\n\n"))
  
    zero_dist <- if(x$dist$zero != "binomial") paste("censored", x$dist$zero, "with log link")
      else paste("binomial with", x$link, "link")
    cat(paste("Zero hurdle model coefficients (", zero_dist, "):\n", sep = ""))
    print.default(format(x$coefficients$zero, digits = digits), print.gap = 2, quote = FALSE)
    if(x$dist$zero == "negbin") cat(paste("Theta =", round(x$theta["zero"], digits), "\n"))
  }
  
  invisible(x)
}

summary.hurdle <- function(object,...)
{
  ## compute z statistics
  kc <- length(object$coefficients$count)
  kz <- length(object$coefficients$zero)
  se <- sqrt(diag(object$vcov))
  coef <- c(object$coefficients$count, object$coefficients$zero)  
  if(object$dist$count == "negbin") {
    coef <- c(coef[1:kc], "Log(theta)" = as.vector(log(object$theta["count"])), coef[(kc+1):(kc+kz)])
    se <- c(se[1:kc], object$SE.logtheta["count"], se[(kc+1):(kc+kz)])
    kc <- kc+1
  }
  if(object$dist$zero == "negbin") {
    coef <- c(coef, "Log(theta)" = as.vector(log(object$theta["zero"])))
    se <- c(se, object$SE.logtheta["zero"])
    kz <- kz+1
  }
  zstat <- coef/se
  pval <- 2*pnorm(-abs(zstat))
  coef <- cbind(coef, se, zstat, pval)
  colnames(coef) <- c("Estimate", "Std. Error", "z value", "Pr(>|z|)")
  object$coefficients$count <- coef[1:kc,,drop = FALSE]
  object$coefficients$zero <- coef[(kc+1):(kc+kz),,drop = FALSE]
  
  ## number of iterations
  object$iterations <- if(!object$separate) tail(na.omit(object$optim$count), 1)
    else tail(na.omit(object$optim$count$count), 1) + tail(na.omit(object$optim$zero$count), 1)
  
  ## delete some slots
  object$residuals <- object$fitted.values <- object$terms <- object$model <- object$y <-
    object$x <- object$levels <- object$contrasts <- object$start <- object$separate <- NULL

  ## return
  class(object) <- "summary.hurdle"
  object
}

print.summary.hurdle <- function(x, digits = max(3, getOption("digits") - 3), ...)
{

  cat("\nCall:", deparse(x$call, width.cutoff = floor(getOption("width") * 0.85)), "\n", sep = "\n")
  
  if(!x$converged) {
    cat("model did not converge\n")
  } else {
    cat(paste("Count model coefficients (truncated ", x$dist$count, " with log link):\n", sep = ""))
    printCoefmat(x$coefficients$count, digits = digits, signif.legend = FALSE)
  
    zero_dist <- if(x$dist$zero != "binomial") paste("censored", x$dist$zero, "with log link")
      else paste("binomial with", x$link, "link")
    cat(paste("Zero hurdle model coefficients (", zero_dist, "):\n", sep = ""))
    printCoefmat(x$coefficients$zero, digits = digits, signif.legend = FALSE)
    
    if(getOption("show.signif.stars") & any(rbind(x$coefficients$count, x$coefficients$zero)[,4] < 0.1))
      cat("---\nSignif. codes: ", "0 '***' 0.001 '**' 0.01 '*' 0.05 '.' 0.1 ' ' 1", "\n")

    if(!is.null(x$theta)) cat(paste("\nTheta:", paste(names(x$theta), round(x$theta, digits), sep = " = ", collapse = ", ")))
    cat(paste("\nNumber of iterations in", x$method, "optimization:", x$iterations, "\n"))
    cat("Log-likelihood:", formatC(x$loglik, digits = digits), "on", x$n - x$df.residual, "Df\n")
  }
  
  invisible(x)
}

terms.hurdle <- function(x, model = c("count", "zero"), ...) {
  x$terms[[match.arg(model)]]
}

model.matrix.hurdle <- function(object, model = c("count", "zero"), ...) {
  model <- match.arg(model)
  if(!is.null(object$x)) rval <- object$x[[model]]
    else if(!is.null(object$model)) rval <- model.matrix(object$terms[[model]], object$model, contrasts = object$contrasts[[model]])
    else stop("not enough information in fitted model to return model.matrix")
  return(rval)
}

predict.hurdle <- function(object, newdata, type = c("response", "prob"), na.action = na.pass, ...)
{
    type <- match.arg(type)

    ## if no new data supplied
    if(missing(newdata)) {
      if(type == "prob") {
        if(!is.null(object$x)) {
	  X <- object$x$count
	  Z <- object$x$zero
	} else if(!is.null(object$model)) {
          X <- model.matrix(object$terms$count, object$model, contrasts = object$contrasts$count)
          Z <- model.matrix(object$terms$zero,  object$model, contrasts = object$contrasts$zero)	
	} else {
	  stop("predicted probabilities cannot be computed with missing newdata")
	}
      } else {
        return(object$fitted.values)
      }
    } else {
      mf <- model.frame(delete.response(object$terms$full), newdata, na.action = na.action, xlev = object$levels)
      X <- model.matrix(delete.response(object$terms$count), mf, contrasts = object$contrasts$count)
      Z <- model.matrix(delete.response(object$terms$zero),  mf, contrasts = object$contrasts$zero)
    }

    phi <- if(object$dist$zero == "binomial") object$linkinv(Z %*% object$coefficients$zero)[,1]
      else exp(Z %*% object$coefficients$zero)[,1]
    p0_zero <- switch(object$dist$zero,
                      "binomial" = 1-phi,
		      "poisson" = exp(-phi),
		      "negbin" = dnbinom(0, size = object$theta["zero"], mu = phi),
		      "geometric" = dnbinom(0, size = 1, mu = phi))

    mu <- exp(X %*% object$coefficients$count)[,1]
    p0_count <- switch(object$dist$count,
		      "poisson" = exp(-mu),
		      "negbin" = dnbinom(0, size = object$theta["count"], mu = mu),
		      "geometric" = dnbinom(0, size = 1, mu = mu))

    
    ## predicted probabilities
    if(type == "prob") {
      if(!is.null(object$y)) y <- object$y
        else if(!is.null(object$model)) y <- model.response(object$model)
	else stop("predicted probabilities cannot be computed for fits with y = FALSE and model = FALSE")

      yUnique <- min(y):max(y)
      nUnique <- length(yUnique)
      rval <- matrix(NA, nrow = length(mu), ncol = nUnique)
      dimnames(rval) <- list(rownames(X), yUnique)
      
      rval[,1] <- p0_zero
      switch(object$dist$count,
             "poisson" = {
               for(i in 2:nUnique) rval[,i] <- (1-p0_zero)/(1-p0_count) * dpois(yUnique[i], lambda = mu)
	     },
	     "negbin" = {
               for(i in 2:nUnique) rval[,i] <- (1-p0_zero)/(1-p0_count) * dnbinom(yUnique[i], mu = mu, size = object$theta["count"])
	     },
	     "geometric" = {
               for(i in 2:nUnique) rval[,i] <- (1-p0_zero)/(1-p0_count) * dnbinom(yUnique[i], mu = mu, size = 1)
	     })
    } else {
      rval <- (1 - p0_zero) * mu / (1 - p0_count)
    }
   
    rval
}

fitted.hurdle <- function(object, ...) {
  object$fitted.values
}

residuals.hurdle <- function(object, type = c("pearson", "response"), ...) {
  type <- match.arg(type)
  if(type == "pearson") return(object$residuals/sqrt(object$fitted.values))
    else return(object$residuals)
}

predprob.hurdle <- function(obj, ...){
    predict(obj, type = "prob", ...)
}

## estfun.hurdle?
