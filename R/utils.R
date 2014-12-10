#' @import multcomp

extrCall <- function (x) 
{
    if (isS4(x)) x@call else x$call
}

extrData <- function (x) 
{
    if (inherits(x, "SuperLearner")) {
        data <- cbind(eval(x$call$X), eval(x$call$Y))
        dimnames(data)[[2]][ncol(data)] <- as.character(x$call$Y[[3]])
        return(data)
    }
    else {
        return(eval(extrCall(x)$data))
    }
}

mgsub <- function (pattern, replacement, x, ...) 
{
    for (i in 1:length(pattern)) x <- gsub(pattern[i], replacement[i], 
        x, ...)
    return(x)
}

neModelEst <- function (formula, family, expData, xFit, ...) 
{
    args <- as.list(match.call())[-1L]
    args$weights <- if (inherits(expData, "weightData")) 
        substitute(attr(expData, "weights"))
    else if (inherits(expData, "impData")) 
        rep(1, nrow(expData))
    if (!is.null(args$xFit)) {
        xFit <- eval(args$xFit)
        class(xFit$formula) <- "Xformula"
        vartypeCheck <- attr(neTerms(xFit$formula), "vartype")
        vartype <- attr(attr(expData, "terms"), "vartype")
        if (!identical(vartype[names(vartypeCheck)], vartypeCheck)) 
            stop("Exposure or covariates of propensity score model and mediator model don't match!")
        dfun <- function(x, fit, data) switch(fit$family$family, 
            gaussian = dnorm(x, mean = predict.glm(fit, newdata = data, 
                type = "response"), sd = sqrt(summary.glm(fit)$dispersion)), 
            binomial = dbinom(as.numeric(x) - 1, size = 1, prob = predict.glm(fit, 
                newdata = data, type = "response")), poisson = dpois(x, 
                lambda = predict.glm(fit, newdata = data, type = "response")))
        denominator <- if (inherits(expData, "weightData")) 
            dfun(expData[, vartype$Xexp[1]], xFit, expData)
        else if (inherits(expData, "impData")) 
            dfun(expData[, vartype$Xexp[2]], xFit, expData)
        args$weights <- eval(args$weights)/denominator
    }
    args$data <- args$expData
    args$expData <- args$xFit <- NULL
    fit <- suppressWarnings(do.call("glm", args))
    mc <- match.call()
    fit$call <- attr(fit, "call") <- mc
    attr(fit, "terms") <- fit$terms
    attr(attr(fit, "terms"), "vartype") <- attr(attr(expData, 
        "terms"), "vartype")
    return(fit)
}
