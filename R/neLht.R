#' Linear hypotheses for natural effect models
#'
#' @description \code{neLht} allows to calculate linear combinations of natural effect model parameter estimates.\cr \code{neEffdecomp} automatically extracts relevant causal parameter estimates from a natural effect model.
#' @param model a fitted natural effect model object.
#' @param ... additional arguments (passed to \code{\link[multcomp]{glht}}).
#' @return An object of class \code{c("neLhtBoot", "neLht", "glht")} (see \code{\link[multcomp]{glht}}). \code{neEffdecomp} returns an object that additionally inherits from class \code{"neEffdecomp"}.
#' 
#' See \code{\link{neLht-methods}} for methods for \code{neLht} objects (and \code{\link[=coef.glht]{glht-methods}} for additional methods for \code{glht} objects).
#' @details \code{neLht} is a wrapper of \code{\link[multcomp]{glht}} and offers the same functionality (see `Details' section of  \code{\link[multcomp]{glht}} for details on argument specification). 
#' It returns objects that inherit from the class \code{"neLht"} in order to make output of their corresponding methods (see \code{\link{neLht-methods}}) more compatible for natural effect models
#' containing bootstrap variance-covariance matrices and standard errors.
#' 
#' \code{neEffdecomp} is a convenience function that automatically extracts causal parameter estimates from a natural effect model
#' and derives natural effect components.
#' That is, natural direct, natural indirect and total causal effect estimates are returned if no exposure-mediator interaction is modelled (i.e. two-way decomposition). 
#' If mediated interaction is allowed for in the natural effect model, there are two ways of decomposing the total effect into (natural) direct and indirect effects components: 
#' either as the sum of the pure direct and the total indirect effect or as the sum of the pure indirect and the total direct effect (i.e. three-way decomposition).
#' In total, five causal effect estimates are returned in this case.
#' @name neLht
#' @note \code{neEffdecomp} is internally called by \code{\link{plot.neModel}} to create confidence interval plots for \code{neModel} objects.
#' @seealso \code{\link{plot.neLht}}, \code{\link{neLht-methods}}, \code{\link[multcomp]{glht}}, \code{\link[=coef.glht]{glht-methods}}, \code{\link{neModel}}, \code{\link{plot.neModel}}
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaffect + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' \donttest{neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                  family = binomial, expData = impData)}\dontshow{neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, family = binomial, expData = impData, nBoot = 2)}
#' 
#' lht <- neLht(neMod, linfct = c("att0 = 0", "att0 + att0:att1 = 0", 
#'                                "att1 = 0", "att1 + att0:att1 = 0", 
#'                                "att0 + att1 + att0:att1 = 0"))
#' summary(lht)
#' 
#' ## or obtain directly via neEffdecomp
#' eff <- neEffdecomp(neMod)
#' summary(eff)
#' @export
NULL

#' Methods for linear hypotheses in natural effect models
#'
#' @description Obtain confidence intervals and statistical tests for linear hypotheses in natural effect models.
#' @param object an object of class \code{neLht}.
#' @param type the type of bootstrap intervals required. The default \code{"norm"} returns normal approximation bootstrap confidence intervals. Currently, only \code{"norm"}, \code{"basic"} and \code{"perc"} are supported (see \code{\link[boot]{boot.ci}}).
#' @param test a function for computing p-values. The default \code{univariate()} does not apply a multiple testing correction. The function \code{adjusted()} allows to correct for multiple testing (see \code{\link[multcomp]{summary.glht}} and \code{\link[multcomp]{adjusted}}) and \code{Chisquare()} allows to test global linear hypotheses.
#' @param ... additional arguments.
#' @inheritParams stats::confint.default
#' @name neLht-methods
#' @details 
#' \code{confint} yields bootstrap confidence intervals. These confidence intervals are internally called via the \code{\link[boot]{boot.ci}} function from the \pkg{boot} package 
#' (and not via the corresponding \code{\link[multcomp]{confint.glht}} function from the \pkg{multcomp} package!).
#' The default confidence level specified in \code{level} (which corresponds to the \code{conf} argument in \code{\link[boot]{boot.ci}}) is 0.95
#' and the default type of bootstrap confidence interval, \code{"norm"}, is based on the normal approximation (for more details see \code{\link[boot]{boot.ci}}).
#'
#' A summary table with large sample tests, similar to that for \code{\link[multcomp]{glht}}, can be obtained using \code{summary}.
#' 
#' In contrast to \code{\link[multcomp]{summary.glht}}, which by default returns \emph{p}-values that are adjusted for multiple testing,
#' the summary function returns unadjusted \emph{p}-values. Adjusted \emph{p}-values can also be obtained by specifying the \code{test} argument 
#' (see \code{\link[multcomp]{adjusted}} for more details).
#' 
#' Global Wald tests considering all linear hypotheses simultaneously (i.e. testing the global null hypothesis) 
#' can be requested by specifying \code{test = Chisqtest()}.
#'
#' See \code{\link[=coef.glht]{glht-methods}} for additional methods for \code{glht} objects.
#'
#' @note \emph{Z}-values in the summary table are simply calculated by dividing the parameter estimate by its corresponding bootstrap standard error. 
#' Corresponding \emph{p}-values in the summary table are only indicative, since the null distribution for each statistic is assumed to be approximately standard normal.
#' Therefore, where possible, it is generally recommend to focus mainly on bootstrap confidence intervals for inference, rather than the provided \emph{p}-values.
#' @seealso \code{\link{neLht}}, \code{\link{plot.neLht}}, \code{\link[multcomp]{glht}}, \code{\link[=coef.glht]{glht-methods}}
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaffect + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' \donttest{neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                  family = binomial, expData = impData)}\dontshow{neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, family = binomial, expData = impData, nBoot = 2)}
#' 
#' lht <- neLht(neMod, linfct = c("att0 = 0", "att0 + att0:att1 = 0", 
#'                                "att1 = 0", "att1 + att0:att1 = 0", 
#'                                "att0 + att1 + att0:att1 = 0"))
#' 
#' ## obtain bootstrap confidence intervals
#' confint(lht)
#' confint(lht, parm = c("att0", "att0 + att0:att1"))
#' confint(lht, parm = 1:2)
#' confint(lht, type = "perc", level = 0.90)
#' 
#' ## summary table
#' summary(lht)
#' 
#' ## summary table with omnibus Chisquare test
#' summary(lht, test = Chisqtest())
NULL

coef.neLhtCI <- function (object, ...) 
{
    attr(object, "coef")
}

#' @rdname neLht-methods
#' @export
confint.neLhtBoot <- function (object, parm, level = 0.95, type = "norm", ...) 
{
    if (missing(parm)) 
        parm <- names(coef(object))
    else if (is.numeric(parm)) 
        parm <- names(coef(object))[parm]
    extractBootci <- function(index, x, level, type) rev(rev(boot::boot.ci(x, 
        conf = level, type = type, t0 = object$linfct[index, 
            ] %*% object$model$bootRes$t0, t = as.vector(object$linfct[index, 
            ] %*% t(object$model$bootRes$t)))[[4]])[1:2])
    ci <- t(sapply(seq.int(nrow(object$linfct)), extractBootci, 
        object$model$bootRes, level = level, type = type))
    dimnames(ci) <- list(rownames(object$linfct), paste0(100 * 
        level, c("% LCL", "% UCL")))
    ci <- ci[parm, ]
    attributes(ci) <- c(attributes(ci), list(level = level, coef = coef(object)[parm], 
        R = attr(object, "R"), type = type))
    class(ci) <- c("neBootCI", "neLhtCI", class(ci))
    return(ci)
}

#' @rdname neLht
#' @export
neEffdecomp <- function (model, ...) 
{
    UseMethod("neEffdecomp")
}

#' @rdname neLht
#' @export
neEffdecomp.neModel <- function (model, ...) 
{
    args <- as.list(match.call())
    orderIA <- function(x) sum(unlist(strsplit(x, "")) == ":") + 
        1
    orderCoef <- sapply(names(coef(model)), orderIA)
    C <- attr(attr(model, "terms"), "vartype")$C
    Xexp <- attr(attr(model, "terms"), "vartype")$Xexp
    ind <- grep(paste(C, collapse = "|"), names(coef(model)))
    orderCoef[ind] <- 1
    K <- diag(length(coef(model)))
    dimnames(K) <- list(names(coef(model)), names(coef(model)))
    ind <- sort(unique(as.vector(sapply(Xexp, grep, dimnames(K)[[1]][orderCoef < 
        2]))))
    K <- K[ind, ]
    if (max(orderCoef) > 1) {
        K <- K[rep(seq.int(nrow(K)), each = 2), ]
        K[, orderCoef > 1] <- rep(c(0, 1), times = 2)
    }
    Klast <- rep(0, length(coef(model)))
    ind <- grep(paste(Xexp, collapse = "|"), names(coef(model)))
    Klast[ind] <- 1
    ind <- grep(paste(C, collapse = "|"), names(coef(model)))
    Klast[ind] <- 0
    K <- rbind(K, Klast)
    extrName <- function(x) names(x[x != 0])
    dimnames(K)[[1]] <- unlist(lapply(apply(K, 1, extrName), 
        function(x) paste(x, collapse = " + ")))
    if (max(orderCoef) == 1) 
        dimnames(K)[[1]] <- c("natural direct effect", "natural indirect effect", 
            "total direct effect")
    else if (max(orderCoef) > 1) 
        dimnames(K)[[1]] <- c("pure direct effect", "total direct effect", 
            "pure indirect effect", "total indirect effect", 
            "total effect")
    effdecomp <- neLht(model, linfct = K)
    attr(effdecomp, "orderCoef") <- orderCoef
    class(effdecomp) <- c("neEffdecomp", class(effdecomp))
    return(effdecomp)
}

#' @rdname neLht
#' @export
neLht <- function (model, ...) 
{
    UseMethod("neLht")
}

#' @rdname neLht
#' @export
neLht.neModelBoot <- function (model, ...) 
{
    lht <- multcomp::glht(model, ...)
    attr(lht, "R") <- model$bootRes$R
    class(lht) <- c("neLhtBoot", "neLht", class(lht))
    return(lht)
}

#' @rdname plot.neLht
#' @export
plot.neEffdecomp <- function (x, level = 0.95, ci.type = "norm", transf = identity, 
    ylabels, yticks.at, ...) 
{
    args <- as.list(match.call())
    if (max(attr(x, "orderCoef")) > 1) 
        if (missing(yticks.at)) 
            args$yticks.at <- c(0, 0.15, 0.5, 0.65, 1)
        else if (missing(yticks.at)) 
            args$yticks.at <- NULL
    args[[1]] <- substitute(plot.neLht)
    eval(as.call(args))
}

#' Confidence interval plots for linear hypotheses in natural effect models
#'
#' @description Confidence interval plots for linear hypotheses in natural effect models.
#' @param x an object of class \code{neLht}.
#' @param ci.type the type of bootstrap intervals required. The default \code{"norm"} returns normal approximation bootstrap confidence intervals. Currently, only \code{"norm"}, \code{"basic"} and \code{"perc"} are supported (see \code{\link[boot]{boot.ci}}).
#' @param transf transformation function to be applied internally on the (linear hypothesis) estimates and their confidence intervals (e.g. \code{exp} for logit or Poisson regression). The default is \code{identity} (i.e. no transformation).
#' @param ylabels character vector containing the labels for the (linear hypothesis) estimates to be plotted on the y-axis.
#' @param yticks.at numeric vector containing the y-coordinates (from 0 to 1) to draw the tick marks for the different estimates and their corresponding confidence intervals.
#' @param ... additional arguments.
#' @inheritParams stats::confint.default
#' @details This function is an adapted version of \code{\link[multcomp]{plot.glht}} from the \pkg{multcomp} package and
#' yields confidence interval plots for each of the linear hypothesis parameters.
#' @seealso \code{\link{neModel}}, \code{\link{neEffdecomp}}, \code{\link{plot.neLht}}, \code{\link{plot.neEffdecomp}}
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaffect + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' \donttest{neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                  family = binomial, expData = impData)}\dontshow{neMod <- neModel(UPB ~ att0 * att1 + gender + educ + age, family = binomial, expData = impData, nBoot = 2)}
#' 
#' lht <- neLht(neMod, linfct = c("att0 = 0", "att0 + att0:att1 = 0", 
#'                                "att1 = 0", "att1 + att0:att1 = 0", 
#'                                "att0 + att1 + att0:att1 = 0"))
#' 
#' ## all pairs return identical output
#' plot(confint(lht), transf = exp)
#' plot(lht, transf = exp)
#' 
#' plot(neEffdecomp(neMod), transf = exp)
#' plot(neMod, transf = exp)
#' 
#' \dontshow{
#'   plot(neEffdecomp(neMod), level = 0.8, transf = exp, ylabels = c("PDE", "TDE", "PIE", "TIE", "TE"), yticks.at = c(0, 0.1, 0.5, 0.6, 1))
#'   plot(neMod, level = 0.8, transf = exp, ylabels = c("PDE", "TDE", "PIE", "TIE", "TE"), yticks.at = c(0, 0.1, 0.5, 0.6, 1))
#'   
#'   lht <- neLht(neMod, linfct = c("att0 = 0"))
#'   summary(lht)
#'   lht <- neLht(neMod, linfct = c("att0 = 0", "att0 + att0:att1 = 2"))
#'   summary(lht)
#' }
#' @export
plot.neLht <- function (x, level = 0.95, ci.type = "norm", transf = identity, 
    ylabels, yticks.at, ...) 
{
    args <- as.list(match.call())[-1L]
    args <- c(substitute(confint), object = substitute(x), type = ci.type, 
        args[names(args) %in% names(formals(confint.neLhtBoot))])
    confint <- eval(as.call(args))
    plot(confint, transf = transf, ylabels = ylabels, yticks.at = yticks.at, 
        ...)
}

plot.neLhtCI <- function (x, transf = identity, ylabels, yticks.at, main, xlab = NULL, 
    ylab = NULL, xlim, ylim, mar, mai, ...) 
{
    args <- as.list(match.call())[-1L]
    if (!grep("confint", args$x) == 1) {
        transf <- eval(args$x[[1]])
        args$x <- args$x[[-1]]
        x <- eval(args$x)
    }
    ci <- transf(x)
    xrange <- c(min(ci[, 1]), max(ci[, 2]))
    yvals <- yticks.at <- if (missing(yticks.at)) 
        nrow(ci):1
    else (nrow(ci) - 1) * (1 - yticks.at) + 1
    if (missing(ylabels)) 
        ylabels <- dimnames(ci)[[1]]
    old.par <- par(no.readonly = TRUE)
    if (all(missing(mai), missing(mar))) {
        mar <- old.par$mar
        mar[2] <- (max(strwidth(ylabels, "inch")) + 0.4) * (old.par$mar/old.par$mai)[1]
    }
    else if (missing(mar)) {
        mar <- mai * (old.par$mar/old.par$mai)[1]
    }
    par(mar = mar)
    null <- transf(0)
    tmp.range <- c(min(null, xrange[1]), max(null, xrange[2]))
    tmp.range <- mean(tmp.range) + c(-1, 1) * 0.6 * diff(tmp.range)
    if (missing(xlab)) 
        xlab <- ""
    if (missing(ylab)) 
        ylab <- ""
    if (missing(xlim)) 
        xlim <- tmp.range
    if (missing(ylim)) 
        ylim <- c(0.5, nrow(ci) + 0.5)
    if (missing(main)) 
        main <- paste0(100 * attr(x, "level"), "% bootstrap CIs")
    plot(c(ci[, 1], ci[, 2]), rep.int(yvals, 2), type = "n", 
        main = main, axes = FALSE, xlab = xlab, ylab = ylab, 
        xlim = xlim, ylim = ylim, ...)
    axis(1, las = 1, ...)
    axis(2, at = yvals, labels = ylabels, las = 1, ...)
    abline(v = null, lty = 2, lwd = 1, ...)
    left <- ci[, 1]
    left[!is.finite(left)] <- min(c(0, xlim[1] * 2))
    right <- ci[, 2]
    right[!is.finite(right)] <- max(c(0, xlim[2] * 2))
    segments(left, yvals, right, yvals, ...)
    points(transf(coef(x)), yvals, pch = 20, ...)
    box()
    par(mar = old.par$mar, mai = old.par$mai)
}

#' @export
print.neLht <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("Linear hypotheses for natural effect models\n---\n")
    coef <- as.matrix(coef(x))
    dimnames(coef)[[2]] <- "Estimate"
    print(coef, digits = digits, ...)
}

#' @export
print.summary.neLht <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    if (inherits(x, "summary.neLhtBoot")) 
        catSE <- "with standard errors based on the non-parametric bootstrap\n"
    else catSE <- NULL
    if (inherits(x$test, "mtest")) {
        cat("Linear hypotheses for natural effect models\n", 
            catSE, "---\n", sep = "")
        if (!identical(x$rhs, rep(0, length(x$rhs)))) 
            dimnames(x$coef.table)[[1]] <- paste(dimnames(x$coef.table)[[1]], 
                "=", x$rhs)
        printCoefmat(x$coef.table, digits = digits, has.Pvalue = TRUE, 
            P.values = TRUE)
        switch(x$test$type, univariate = cat("(Univariate p-values reported)"), 
            `single-step` = cat("(Adjusted p-values reported -- single-step method)"), 
            Shaffer = cat("(Adjusted p-values reported -- Shaffer method)"), 
            Westfall = cat("(Adjusted p-values reported -- Westfall method)"), 
            cat("(Adjusted p-values reported --", x$test$type, 
                "method)"))
        cat("\n\n")
    }
    else if (inherits(x$test, "gtest")) {
        cat("Global linear hypothesis test for natural effect models\n", 
            catSE, "---\n", sep = "")
        pr <- data.frame(x$test$SSH, x$test$df[1], x$test$pval)
        names(pr) <- c("Chisq", "DF", "Pr(>Chisq)")
        print(pr, digits = digits, ...)
    }
}

#' @rdname neLht-methods
#' @export
summary.neLht <- function (object, test = univariate(), ...) 
{
    class(object) <- "glht"
    summary <- summary(object, test = test, ...)
    pq <- summary$test
    if (inherits(pq, "mtest")) {
        coef.table <- cbind(pq$coefficients, pq$sigma, pq$tstat, 
            pq$pvalues)
        dimnames(coef.table) <- list(names(coef(object)), c("Estimate", 
            "Std. Error", "z value", "Pr(>|z|)"))
        summary$coef.table <- coef.table
    }
    class(summary) <- c("summary.neLhtBoot", "summary.neLht", 
        class(summary))
    return(summary)
}
