#' Methods for natural effect models
#'
#' @description Confidence intervals and statistical tests for natural effect models.
#' @param object a fitted natural effect model object.
#' @param ... additional arguments (see \code{\link[boot]{boot.ci}} for \code{confint} or \code{\link[stats]{summary.glm}} for \code{summary}).
#' @inheritParams stats::confint.default
#' @inheritParams neLht-methods
#' @details \code{vcov} returns the variance covariance matrix calculated from the bootstrap samples stored in\cr \code{object$bootRes} (see \code{\link{neModel}}).
#'
#' \code{confint} yields bootstrap confidence intervals. These confidence intervals are internally called via the \code{\link[boot]{boot.ci}} function from the \pkg{boot} package.
#' The default confidence level specified in \code{level} (which corresponds to the \code{conf} argument in \code{\link[boot]{boot.ci}}) is 0.95
#' and the default type of bootstrap confidence interval, \code{"norm"}, is based on the normal approximation (for more details see \code{\link[boot]{boot.ci}}).
#'
#' A summary table with large sample tests, similar to that for \code{\link[stats]{glm}} output, can be obtained using \code{summary}.
#'
#' @name neModel-methods
#' @note \emph{Z}-values in the summary table are simply calculated by dividing the parameter estimate by its corresponding bootstrap standard error. 
#' Corresponding \emph{p}-values in the summary table are only indicative, since the null distribution for each statistic is assumed to be approximately standard normal.
#' Therefore, where possible, it is generally recommend to focus mainly on bootstrap confidence intervals for inference, rather than the provided \emph{p}-values.
#' @seealso \code{\link{neModel}}, \code{\link{plot.neModel}}
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaffect + educ + gender + age, 
#'                     family = binomial, data = UPBdata)
#' \donttest{neMod <- neModel(UPB ~ att0 * att1 + educ + gender + age, 
#'                  family = binomial, expData = impData)}\dontshow{neMod <- neModel(UPB ~ att0 * att1 + educ + gender + age, family = binomial, expData = impData, nBoot = 2)}
#' 
#' ## extract coefficients
#' coef(neMod)
#' 
#' ## extract (bootstrap) variance-covariance matrix
#' vcov(neMod)
#' 
#' ## obtain bootstrap confidence intervals
#' confint(neMod)
#' confint(neMod, parm = c("att0"))
#' confint(neMod, type = "perc", level = 0.90)
#' 
#' ## summary table
#' summary(neMod)
NULL

#' @rdname neModel-methods
#' @export
coef.neModel <- function (object, ...) 
{
    coef(object$neModelFit)
}

#' @rdname neModel-methods
#' @export
confint.neModelBoot <- function (object, parm, level = 0.95, type = "norm", ...) 
{
    if (missing(parm)) 
        parm <- names(coef(object))
    else if (is.numeric(parm)) 
        parm <- names(coef(object))[parm]
    extractBootci <- function(index, x, level, type) rev(rev(boot::boot.ci(x, 
        conf = level, type = type, index = index)[[4]])[1:2])
    ci <- t(sapply(1:length(object$bootRes$t0), extractBootci, 
        object$bootRes, level = level, type = type))
    dimnames(ci) <- list(names(coef(object)), paste0(100 * level, 
        c("% LCL", "% UCL")))
    ci <- ci[parm, ]
    attributes(ci) <- c(attributes(ci), list(level = level, coef = coef(object)[parm], 
        R = object$bootRes$R, type = type))
    class(ci) <- c("neBootCI", "neModelBootCI", class(ci))
    return(ci)
}

df.residual.neModel <- function (object, ...) 
{
    df.residual(object$neModelFit, ...)
}

model.matrix.neModel <- function (object, ...) 
{
    model.matrix(object$neModelFit, ...)
}

#' Fit a natural effect model
#'
#' @description \code{neModel} is used to fit a natural effect model on the expanded dataset.
#' @param formula a \code{\link[stats]{formula}} object providing a symbolic description of the natural effect model.
#' @param expData the expanded dataset (of class \code{"\link{expData}"}).
#' @param xFit fitted model object representing a model for the exposure (used for inverse treatment probability weighting).
#' @param nBoot number of bootstrap replicates (see \code{R} argument of \code{\link[boot]{boot}}).
#' @param ncpus integer: number of processes to be used in parallel operation: typically one would chose this to the number of available CPUs (see details).
#' @param progress logical value indicating whether or not a progress bar should be displayed. Progress bars are automatically disabled for multicore processing.
#' @param ... additional arguments (passed to \code{\link[stats]{glm}}).
#' @inheritParams stats::glm
#' @inheritParams boot::boot
#' @return An object of class \code{"\link[=neModel-methods]{neModel}"} consisting of a list of 2 objects:
#' \item{\code{neModelFit}}{the fitted natural model object (of class \code{"\link[stats]{glm}"}) with downwardly biased standard errors}
#' \item{\code{bootRes}}{the bootstrap results (of class \code{"\link[boot]{boot}"})}
#' See \code{\link{neModel-methods}} for methods for \code{neModel} objects.
#' @details This function is a wrapper for \code{\link[stats]{glm}}, providing unbiased bootstrap standard errors for the parameter estimates.
#'
#' The \code{formula} argument requires to be specified in function of the variables from the expanded dataset (specified in \code{expData}) whose corresponding parameters index the direct and indirect effect.
#' Stratum-specific natural effects can be estimated by additionally modeling the relation between the outcome and baseline covariates.
#' If the set of baseline covariates adjusted for in the \code{formula} argument is not sufficient to control for confounding (e.g. when fitting a population-average natural effect model),
#' an adequate model for the exposure (conditioning on a sufficient set of baseline covariates) should be specified in the \code{xFit} argument.
#' In this case, such a model for the exposure distribution is needed to weight by the reciprocal of the probability (density) of the exposure (i.e. inverse probability weighting) in order to adjust for confounding.
#' Just as for ratio-of-mediator probability weighting (see paragraph below), this kind of weighting is done internally.
#'
#' In contrast to \code{\link[stats]{glm}}, the \code{expData} argument (rather than \code{data} argument) requires specification of a data frame that inherits from the class \code{"\link{expData}"},
#' which contains additional information about e.g. the fitted working model, the variable types or terms of this working model 
#' and possibly ratio-of-mediator probability weights.
#' The latter are automatically extracted from the \code{\link{expData}} object and weighting is done internally.
#'
#' As the default \code{\link[stats]{glm}} standard errors fail to reflect the uncertainty inherent to the working model(s) (i.e. either a model for the mediator or an imputation model for the outcome and possibly a model for the exposure),
#' bootstrap standard errors are calculated (using the \code{\link[boot]{boot}} function from the \pkg{boot} package). The bootstrap procedure entails refitting these working models on each bootstrap sample, reconstructing the expanded dataset and subsequently refitting the specified natural effect model on this dataset.
#' In order to obtain stable standard errors, the number of bootstrap samples (specified via the \code{nBoot} argument) should be chosen relatively high. Although the default is 1000, it is recommended to use at least 500 bootstrap samples.
#' 
#' To speed up the bootstrap procedure, parallel processing can be used by specifying the desired type of parallel operation via the \code{parallel} argument (for more details, see \code{\link[boot]{boot}}).
#' The number of parallel processes (\code{ncpus}) is suggested to be specified explicitly (its default is 1, unless the global option \code{options("boot.cpus")} is specified). 
#' The function \code{\link[parallel]{detectCores}} from the \pkg{parallel} package can be helpful at determining the number of available cores (although this may not always correspond to the number of \emph{allowed} cores).
#'
#' @note It is important to note that the original mediator(s) should not be specified in the \code{formula} argument, as the natural indirect effect in natural effect models
#' should by captured solely by parameter(s) corresponding to the auxiliary hypothetical variable \emph{x*} in the expanded dataset (see \code{\link{expData}}).
#'
#' @seealso \code{\link{neModel-methods}}, \code{\link{plot.neModel}}, \code{\link{neImpute}}, \code{\link{neWeight}}, \code{\link{neLht}}, \code{\link{neEffdecomp}}
#'
#' @examples
#' data(UPBdata)
#' 
#' ## weighting-based approach
#' weightData <- neWeight(negaffect ~ att + gender + educ + age, 
#'                        data = UPBdata)
#' 
#' # stratum-specific natural effects
#' \donttest{weightFit1 <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                       family = binomial, expData = weightData)}\dontshow{weightFit1 <- neModel(UPB ~ att0 * att1 + gender + educ + age, family = binomial, expData = weightData, nBoot = 2)}
#' summary(weightFit1)
#' 
#' # population-average natural effects
#' expFit <- glm(att ~ gender + educ + age, data = UPBdata)
#' \donttest{weightFit2 <- neModel(UPB ~ att0 * att1, family = binomial, 
#'                       expData = weightData, xFit = expFit)}\dontshow{weightFit2 <- neModel(UPB ~ att0 * att1, family = binomial, expData = weightData, xFit = expFit, nBoot = 2)}
#' summary(weightFit2)
#' 
#' ## imputation-based approach
#' impData <- neImpute(UPB ~ att * negaffect + gender + educ + age, 
#'                     family = binomial, data = UPBdata)
#' 
#' # stratum-specific natural effects
#' \donttest{impFit1 <- neModel(UPB ~ att0 * att1 + gender + educ + age, 
#'                    family = binomial, expData = impData)
#' summary(impFit1)}
#' 
#' # population-average natural effects
#' \donttest{impFit2 <- neModel(UPB ~ att0 * att1, family = binomial, 
#'                    expData = impData, xFit = expFit)}\dontshow{impFit2 <- neModel(UPB ~ att0 * att1, family = binomial, expData = impData, xFit = expFit, nBoot = 2)}
#' summary(impFit2)
#' 
#' \dontshow{# check with vgam (VGAM package)
#' library(VGAM)
#' weightData <- neWeight(negaffect ~ att + gender + educ + age, family = "gaussianff", data = UPBdata, FUN = vgam)
#' impData <- neImpute(UPB ~ att + negaffect + gender + educ + age, family = "binomialff", data = UPBdata, FUN = vgam)
#' # debug(neModel)
#' weightFit <- neModel(UPB ~ att0 + att1 + gender + educ + age, family = binomial, expData = weightData, nBoot = 2)
#' impFit <- neModel(UPB ~ att0 + att1 + gender + educ + age, family = binomial, expData = impData, nBoot = 2)
#' summary(weightFit)
#' summary(impFit)
#'
#' # warning!
#' impFit <- neModel(UPB ~ att0 * att1 + gender + educ + age, family = binomial, expData = impData, nBoot = 2)
#'
#' # inverse propensity score weighting
#' expFit <- glm(att ~ gender + educ + age, data = UPBdata)
#' impFit <- neModel(UPB ~ att0 + att1, family = binomial, expData = impData, xFit = expFit, nBoot = 2)}
#' @references
#' Lange, T., Vansteelandt, S., & Bekaert, M. (2012). A Simple Unified Approach for Estimating Natural Direct and Indirect Effects. \emph{American Journal of Epidemiology}, \bold{176}(3), 190-195.
#' 
#' Vansteelandt, S., Bekaert, M., & Lange, T. (2012). Imputation Strategies for the Estimation of Natural Direct and Indirect Effects. \emph{Epidemiologic Methods}, \bold{1}(1), Article 7.
#' 
#' Loeys, T., Moerkerke, B., De Smet, O., Buysse, A., Steen, J., & Vansteelandt, S. (2013). Flexible Mediation Analysis in the Presence of Nonlinear Relations: Beyond the Mediation Formula. \emph{Multivariate Behavioral Research}, \bold{48}(6), 871-894.
#' @export
neModel <- function (formula, family = gaussian, expData, xFit, nBoot = 1000, 
    parallel = c("no", "multicore", "snow"), ncpus = getOption("boot.ncpus", 
        1L), progress = TRUE, ...) 
{
    args <- as.list(match.call())
    args[[1]] <- substitute(neModelEst)
    args <- c(args[[1]], args[names(args) %in% names(formals(neModelEst))])
    neModelFit <- eval(as.call(args))
    fit <- attr(expData, "model")
    parallel <- match.arg(parallel)
    if (parallel != "no") {
        if (progress) {
            message("Progress bar disabled for parallel processing.")
            progress <- FALSE
        }
    }
    if (progress) {
        env <- environment()
        counter <- 0
        progbar <- txtProgressBar(min = 0, max = nBoot, style = 3)
    }
    if (sum(sapply(attr(terms(expData), "vartype")$M, grepl, 
        labels(terms(neModelFit))))) 
        stop("The original mediator variables should not be included in the natural effect model!")
    if (inherits(expData, "impData")) {
        formulaImpData <- extrCall(attr(expData, "model"))$formula
        if (is.null(formulaImpData)) 
            formulaImpData <- attr(expData, "call")$formula
        termsImpData <- mgsub(dimnames(attr(terms(expData), "factors"))[[1]], 
            all.vars(formulaImpData), labels(terms(expData)), 
            fixed = TRUE)
        termsImpData <- mgsub(unlist(attr(terms(expData), "vartype")[c("X", 
            "M")]), attr(terms(expData), "vartype")$Xexp, termsImpData, 
            fixed = TRUE)
        termsNeModel <- mgsub(dimnames(attr(terms(neModelFit), 
            "factors"))[[1]], all.vars(neModelFit$formula), labels(terms(neModelFit)), 
            fixed = TRUE)
        if (sum(!termsNeModel %in% termsImpData)) 
            warning("The imputation model should at least reflect the structure of the natural effect model! This may cause some of the effect estimates to be attenuated. Please consult the terms of the models 'x' using 'labels(terms(x))'.")
    }
    neModelBootfun <- function(data, ind, fit, expData, neModelFit) {
        bootData <- data[ind, ]
        FUN <- extrCall(fit)[[1]]
        call1 <- as.list(match.call(eval(FUN), extrCall(fit)))
        if (inherits(fit, "SuperLearner")) 
            call1$X[[2]] <- call1$Y[[2]] <- substitute(bootData)
        else call1$data <- substitute(bootData)
        bootExpFit <- eval(as.call(call1))
        nRep <- nrow(expData)/nrow(data)
        indExp <- as.vector(sapply(ind, function(x) (x * nRep - 
            (nRep - 1)):(x * nRep)))
        bootExpData <- expData[indExp, ]
        bootExpData <- eval(attr(expData, "call")[[1]])(bootExpFit, 
            skipExpand = TRUE, expData = bootExpData)
        if (!is.null(extrCall(neModelFit)$xFit)) {
            call3 <- as.list(eval(extrCall(neModelFit)$xFit)$call)
            call3$data <- substitute(bootData)
            bootxFit <- eval(as.call(call3))
        }
        call4 <- as.list(attr(neModelFit, "call"))
        call4$expData <- substitute(bootExpData)
        if (!is.null(extrCall(neModelFit)$xFit)) 
            call4$xFit <- substitute(bootxFit)
        bootFit <- eval(as.call(call4))
        if (progress) {
            curVal <- get("counter", envir = env)
            assign("counter", curVal + 1, envir = env)
            setTxtProgressBar(get("progbar", envir = env), curVal + 
                1)
        }
        return(coef(bootFit))
    }
    bootRes <- boot::boot(data = extrData(fit), statistic = neModelBootfun, 
        R = nBoot, fit = fit, expData = expData, neModelFit = neModelFit, 
        parallel = parallel, ncpus = ncpus)
    neModel <- list(neModelFit = neModelFit, bootRes = bootRes)
    attr(neModel, "terms") <- neModelFit$terms
    attr(attr(neModel, "terms"), "vartype") <- attr(attr(expData, 
        "terms"), "vartype")
    class(neModel) <- c("neModelBoot", "neModel")
    return(neModel)
}

#' Confidence interval plots for natural effect estimates
#'
#' @description Obtain effect decomposition confidence interval plots for natural effect models.
#' @param x a fitted natural effect model object.
#' @inheritParams plot.neLht
#' @details This function yields confidence interval plots for the natural effect parameters of interest.
#' These causal parameter estimates are first internally extracted from the \code{neModel} object by applying\cr\code{\link{neEffdecomp}}.
#'
#' @examples
#' data(UPBdata)
#' 
#' impData <- neImpute(UPB ~ att * negaffect + educ + gender + age, 
#'                     family = binomial, data = UPBdata)
#' \donttest{neMod <- neModel(UPB ~ att0 * att1 + educ + gender + age, 
#'                  family = binomial, expData = impData)}\dontshow{neMod <- neModel(UPB ~ att0 * att1 + educ + gender + age, family = binomial, expData = impData, nBoot = 2)}
#'
#' plot(neMod)
#' plot(neMod, transf = exp, 
#'      ylabels = c("PDE", "TDE", "PIE", "TIE", "TE"))
#' @export
plot.neModel <- function (x, level = 0.95, ci.type = "norm", transf = identity, 
    ylabels, yticks.at, ...) 
{
    args <- as.list(match.call())
    args[[1]] <- substitute(plot)
    args$x <- substitute(neEffdecomp(x))
    eval(as.call(args))
}

#' @export
print.neBootCI <- function (x, ...) 
{
    attr <- attributes(x)
    attributes(x)[c("level", "coef", "R", "type", "class")] <- NULL
    print.default(x)
    cat("---\n")
    cat(switch(attr$type, norm = "normal approximation bootstrap", 
        basic = "basic bootstrap", stud = "studentized bootstrap", 
        perc = "bootstrap percentile", bca = "adjusted bootstrap percentile (BCa)"), 
        " confidence intervals\nbased on ", attr$R, " bootstrap samples", 
        sep = "")
}

#' @export
print.summary.neModel <- function (x, digits = max(3, getOption("digits") - 3), ...) 
{
    cat("Natural effect model\n")
    cat("with standard errors based on the non-parametric bootstrap\n---\n")
    cat("Exposure:", x$terms$X, "\nMediator(s):", paste(x$terms$M, 
        collapse = ", "), "\n---\n")
    cat("Parameter estimates:\n")
    printCoefmat(x$coef.table, digits = digits, has.Pvalue = TRUE, 
        P.values = TRUE)
}

#' @rdname neModel-methods
#' @export
summary.neModel <- function (object, ...) 
{
    summary <- summary(object$neModelFit)
    se <- sqrt(diag(vcov(object)))
    zvalue <- coef(object)/se
    pvalue <- 2 * pnorm(-abs(zvalue))
    coef.table <- cbind(coef(object), se, zvalue, pvalue)
    dimnames(coef.table) <- list(names(coef(object)), c("Estimate", 
        "Std. Error", "z value", "Pr(>|z|)"))
    summary$coef.table <- coef.table
    summary$terms <- attr(attr(object, "terms"), "vartype")
    class(summary) <- "summary.neModel"
    return(summary)
}

#' @rdname neModel-methods
#' @export
vcov.neModelBoot <- function (object, ...) 
{
    covmat <- cov(object$bootRes$t)
    dimnames(covmat) <- dimnames(summary.glm(object$neModelFit)$cov.scaled)
    return(covmat)
}
