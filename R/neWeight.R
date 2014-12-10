#' Expand the dataset and calculate ratio-of-mediator probability weights
#'
#' @description This function both expands the data along hypothetical exposure values and calculates ratio-of-mediator probability weights.
#' @param object an object used to select a method.
#' @param ... additional arguments.
#' @return A data frame of class \code{c("data.frame", "expData", "weightData"))}. See \code{\link{expData}} for its structure.
#' @details Generic function that both expands the data along hypothetical exposure values and
#' calculates ratio-of-mediator probability weights 
#' 
#' \deqn{\frac{\hat P(M_i \vert X_i = x^*, C_i)}{\hat P(M_i \vert X_i = x, C_i)}}
#' 
#' for each observation unit \emph{i} in this expanded dataset in a single run.
#' These weights are ratios of probabilities or probability densities from the mediator model distribution, which can be specified either externally as a fitted model object (\code{\link{neWeight.default}})
#' or internally (\code{\link{neWeight.formula}}).
#' @seealso \code{\link{neWeight.default}}, \code{\link{neWeight.formula}}, \code{\link{expData}}
#' @references
#' Lange, T., Vansteelandt, S., & Bekaert, M. (2012). A Simple Unified Approach for Estimating Natural Direct and Indirect Effects. \emph{American Journal of Epidemiology}, \bold{176}(3), 190-195.
#' @export
neWeight <- function (object, ...) 
{
    UseMethod("neWeight")
}

#' Expand the dataset and calculate ratio-of-mediator probability weights
#'
#' @description This function both expands the data along hypothetical exposure values and calculates ratio-of-mediator probability weights.
#' @param object fitted model object representing the mediator model.
#' @param formula a \code{\link[stats]{formula}} object providing a symbolic description of the mediator model. Redundant if already specified in call for fitted model specified in \code{object} (see details).
#' @inheritParams neImpute.default
#' @return A data frame of class \code{c("data.frame", "expData", "weightData"))}. See \code{\link{expData}} for its structure.
#' @details The calculated weights are ratios of fitted probabilities or probability densities from the distribution of the mediator model.
#' This model needs to be specified as a fitted object in the \code{object} argument.
#'
#' If the model-fitting function used to fit the mediator model does not require specification of a \code{formula} or \code{data} argument,
#' these need to be specified explicitly in order to enable \code{neWeight.default} to extract pointers to variable types relevant for mediation analysis.
# (also see \code{\link{neTerms}}).
#'
#' Whether a \code{\link[stats]{formula}} is specified externally (in the call for the fitted mediator model object which is specified in \code{object}) or internally (via the \code{formula} argument),
#' it always needs to be of the form \code{M ~ X + C1 + C2}, with predictor variables entered in the following prespecified order:
#' \enumerate{
#'  \item exposure \code{X}: The first predictor is coded as exposure or treatment.
#'  \item baseline covariates \code{C}: All remaining predictor variables are automatically coded as baseline covariates.
#' }
#'
#' In contrast to imputation models with categorical exposures, additional arguments need to be specified if the exposure is continuous.
#' All of these additional arguments are related to the sampling procedure for the exposure.
#'
#' Whereas the number of replications \code{nRep} for categorical variables equals the number of levels for the exposure coded as a factor (i.e. the number of hypothetical exposure values), the number of desired replications needs to be specified explicitly for continuous exposures.
#' Its default is 5.
#'
#' If \code{xFit} is left unspecified, the hypothetical exposure levels are automatically sampled from a linear model for the exposure, conditional on a linear combination of all covariates.
#' If one wishes to use another model for the exposure, this default model specification can be overruled by referring to a fitted model object in the \code{xFit} argument.
#' Misspecification of this sampling model does not induce bias in the estimated coefficients and standard errors of the natural effect model.
#'
#' The \code{xSampling} argument allows to specify how the hypothetical exposure levels should be sampled from the conditional exposure distribution (which is either entered explicitly using the \code{xFit} argument or fitted automatically as described in the previous paragraph).
#' The \code{"random"} option randomly samples \code{nRep} draws from the exposure distribution, whereas the \code{"quantiles"} option (default) samples \code{nRep} quantiles at equal-sized probability intervals. Only the latter hence yields fixed exposure levels given \code{nRep} and \code{xFit}. \cr\cr
#' In order to guarantee that the entire support of the distribution is being sampled (which might be a concern if \code{nRep} is chosen to be small), the default lower and upper sampled quantiles are the 5th and 95th percentiles.
#' The intermittent quantiles correspond to equal-sized probability intervals. So, for instance, if \code{nRep = 4}, then the sampled quantiles will correspond to probabilities 0.05, 0.35, 0.65 and 0.95.
#' These default 'outer' quantiles can be changed by specifying the \code{percLim} argument accordingly. By specifying \code{percLim = NULL}, the standard quantiles will be sampled (e.g., 0.2, 0.4, 0.6 and 0.8 if \code{nRep = 4}).
#'
#' @seealso \code{\link{neWeight}}, \code{\link{neWeight.formula}}, \code{\link{expData}}
#' @examples
#' data(UPBdata)
#' 
#' ## example using glm
#' fit.glm <- glm(negaffect ~ att + gender + educ + age, data = UPBdata)
#' weightData <- neWeight(fit.glm, nRep = 2)
#' weights <- attr(weightData, "weights")
#' head(cbind(weightData, weights))
#' 
#' ## example using vglm (yielding identical results as with glm)
#' library(VGAM)
#' fit.vglm <- vglm(negaffect ~ att + gender + educ + age, 
#'                  family = gaussianff, data = UPBdata)
#' weightData2 <- neWeight(fit.vglm, nRep = 2)
#' weights2 <- attr(weightData2, "weights")
#' head(cbind(weightData2, weights2))
#'
#' \dontshow{fit1 <- glm(negaffect ~ att + gender + educ + age, data = UPBdata)
#' expData1 <- neWeight(fit1)
#' w1 <- attr(expData1, "weights")
#' expData1f <- neWeight(negaffect ~ att + gender + educ + age, data = UPBdata)
#' w1f <- attr(expData1f, "weights")
#' head(expData1); head(expData1f)
#' head(w1); head(w1f)
#'
#' # test vglm (vglm is also vgam class, but not other way around!)
#' library(VGAM)
#' fit1b <- vgam(negaffect ~ att + gender + educ + age, family = gaussianff, data = UPBdata)
#' expData1b <- neWeight(fit1b)
#' head(attr(expData1, "weights")); head(attr(expData1b, "weights"))
#' fit1b <- vgam(negaffect ~ s(att) + gender + educ + age, family = gaussianff, data = UPBdata)
#' expData1b <- neWeight(fit1b)
#' head(attr(expData1, "weights")); head(attr(expData1b, "weights"))
#' expData1bf <- neWeight(negaffect ~ s(att) + gender + educ + age, FUN = vgam, family = gaussianff, data = UPBdata)
#' head(attr(expData1b, "weights")); head(attr(expData1bf, "weights"))
#' ##
#'
#' UPBdata$negaffect2 <- cut(UPBdata$negaffect, breaks = 2, labels = c("low", "high"))
#' fit2 <- glm(negaffect2 ~ att + gender + educ + age, family = binomial, data = UPBdata)
#' expData2 <- neWeight(fit2)
#' w2 <- attr(expData2, "weights")
#' expData2f <- neWeight(negaffect2 ~ att + gender + educ + age, family = binomial, data = UPBdata)
#' w2f <- attr(expData2f, "weights")
#' head(expData2); head(expData2f)
#' head(w2); head(w2f)
#'
#' # test vglm
#' fit2b <- vgam(negaffect2 ~ att + gender + educ + age, family = binomialff, data = UPBdata)
#' expData2b <- neWeight(fit2b)
#' head(attr(expData2, "weights")); head(attr(expData2b, "weights"))
#' fit2b <- vgam(negaffect2 ~ s(att) + gender + educ + age, family = binomialff, data = UPBdata)
#' expData2b <- neWeight(fit2b)
#' head(attr(expData2, "weights")); head(attr(expData2b, "weights"))
#' expData2bf <- neWeight(negaffect2 ~ s(att) + gender + educ + age, FUN = vgam, family = binomialff, data = UPBdata)
#' head(attr(expData2b, "weights")); head(attr(expData2bf, "weights"))
#' ##
#'
#' UPBdata$negaffect3 <- cut(UPBdata$negaffect, breaks = 3, labels = c("low", "moderate", "high"))
#' UPBdata$negaffect3 <- as.numeric(UPBdata$negaffect3)
#' fit3 <- glm(negaffect3 ~ att + gender + educ + age, family = "poisson", data = UPBdata)
#' expData3 <- neWeight(fit3)
#' w3 <- attr(expData3, "weights")
#' expData3f <- neWeight(negaffect3 ~ att + gender + educ + age, family = poisson, data = UPBdata)
#' w3f <- attr(expData3f, "weights")
#' head(expData3); head(expData3f)
#' head(w3); head(w3f)
#'
#' # test vglm
#' fit3b <- vgam(negaffect3 ~ att + gender + educ + age, family = poissonff, data = UPBdata)
#' expData3b <- neWeight(fit3b)
#' head(attr(expData3, "weights")); head(attr(expData3b, "weights"))
#' fit3b <- vgam(negaffect3 ~ s(att) + gender + educ + age, family = poissonff, data = UPBdata)
#' expData3b <- neWeight(fit3b)
#' head(attr(expData3, "weights")); head(attr(expData3b, "weights"))
#' expData3bf <- neWeight(negaffect3 ~ s(att) + gender + educ + age, FUN = vgam, family = poissonff, data = UPBdata)
#' head(attr(expData3b, "weights")); head(attr(expData3bf, "weights"))}
#' @export
neWeight.default <- function (object, formula, data, nRep = 5, xSampling = c("quantiles", 
    "random"), xFit, percLim = c(0.05, 0.95), ...) 
{
    args <- as.list(match.call())[-1L]
    nMed <- args$nMed <- 1
    fit <- object
    if (!isTRUE(args$skipExpand)) {
        args$data <- if (missing(data)) {
            extrCall(fit)$data
        }
        else substitute(data)
        if (missing(formula)) 
            formula <- extrCall(fit)$formula
        class(formula) <- c("Mformula", "formula")
        vartype <- args$vartype <- attr(neTerms(formula, Y = NA, 
            nMed = nMed, joint = joint), "vartype")
        xFact <- try(is.factor(model.frame(formula, eval(args$data))[, 
            attr(vartype, "xasis")]), silent = TRUE)
        if (!is.logical(xFact)) 
            xFact <- is.factor(model.frame(fit)[, attr(vartype, 
                "xasis")])
        if (xFact) {
            dontapply <- c("nRep", "xSampling", "xFit", "percLim")
            ind <- !sapply(args[dontapply], is.null)
            ind <- ind[which(ind)]
            if (length(ind) > 1) 
                warning(gettextf("The arguments %s don't apply when the exposure is categorical!", 
                  paste(paste0("'", names(ind), "'"), collapse = ", ")))
            else if (length(ind) == 1) 
                warning(gettextf("The argument '%s' doesn't apply when the exposure is categorical!", 
                  names(ind)))
        }
        args$xSampling <- xSampling <- match.arg(xSampling)
        if (all(!missing(percLim), xSampling == "random")) 
            warning("The percentile limits that have been specified in 'percLim' do not apply when 'xSampling' is specified as 'random'!")
        joint <- args$joint <- TRUE
        if (missing(nRep)) 
            args$nRep <- substitute(nRep)
        if (missing(percLim)) 
            args$percLim <- percLim
        expData <- do.call("expandData", c(x = substitute(vartype$X), 
            args))
        nExp <- ifelse(joint, 1, nMed)
        tmp <- vartype$Xexp <- paste0(vartype$X, seq(0, nExp))
        names(expData)[sapply(paste0("aux", seq(0, nExp)), grep, 
            names(expData))] <- tmp
        other <- names(expData)[!is.element(names(expData), c("id", 
            vartype$X, vartype$Xexp))]
        expData <- expData[, c("id", vartype$Xexp, other)]
    }
    else {
        if (is.null(args$expData)) 
            stop("expData is missing!")
        expData <- eval.parent(args$expData)
        attr <- attributes(expData)
        vartype <- attr(attr$terms, "vartype")
    }
    family <- if (is.null(extrCall(fit)$family)) 
        formals(eval(extrCall(fit)[[1]]))$family
    else extrCall(fit)$family
    family <- c("gaussian", "binomial", "poisson")[mapply(function(x, 
        y) grepl(y, x), as.character(family), c("gaussian", "binomial", 
        "poisson"))]
    dispersion <- if (inherits(fit, "vglm")) 
        fit@misc$dispersion
    else summary(fit)$dispersion
    dfun <- function(x) switch(family, gaussian = dnorm(x, mean = predict(fit, 
        newdata = expData, type = "response"), sd = sqrt(dispersion)), 
        binomial = dbinom(as.numeric(x) - 1, size = 1, prob = predict(fit, 
            newdata = expData, type = "response")), poisson = dpois(x, 
            lambda = predict(fit, newdata = expData, type = "response")))
    expData[, vartype$X] <- expData[, vartype$Xexp[2]]
    weightsNum <- dfun(expData[, vartype$M])
    expData[, vartype$X] <- expData[, vartype$Xexp[1]]
    weightsDenom <- dfun(expData[, vartype$M])
    expData <- expData[, -ncol(expData)]
    if (!isTRUE(args$skipExpand)) {
        attributes(expData) <- c(attributes(expData), list(model = fit, 
            call = match.call(), terms = neTerms(formula, nMed, 
                joint), weights = weightsNum/weightsDenom))
        attr(attr(expData, "terms"), "vartype") <- vartype
        class(expData) <- c(class(expData), "expData", "weightData")
    }
    else {
        attributes(expData) <- attr
        attr(expData, "model") <- fit
        attr(expData, "weights") <- weightsNum/weightsDenom
    }
    return(expData)
}

#' Expand the dataset and calculate ratio-of-mediator probability weights
#'
#' @description This function both expands the data along hypothetical exposure values and calculates ratio-of-mediator probability weights.
#' @param object a \code{\link[stats]{formula}} object providing a symbolic description of the mediator model (see details).
#' @inheritParams neImpute.formula
#' @return A data frame of class \code{c("data.frame", "expData", "weightData"))}. See \code{\link{expData}} for its structure.
#' @details The calculated weights are ratios of fitted probabilities or probability densities from the distribution of the mediator model.
#' This model is fitted internally by extracting information from the arguments \code{object}, \code{family}, \code{data}, \code{FUN} and \code{...}.
#'
#' For mediation model specification via the \code{object} argument, use a \code{\link[stats]{formula}} of the form\cr\code{M ~ X + C1 + C2},
#' with predictor variables entered in the following prespecified order:
#' \enumerate{
#'  \item exposure \code{X}: The first predictor is coded as exposure or treatment.
#'  \item baseline covariates \code{C}: All remaining predictor variables are automatically coded as baseline covariates.
#' }
#'
#' The type of mediator model can be defined by specifying an appropriate model-fitting function via the \code{FUN} argument (its default is \code{\link[stats]{glm}}).
#' This method can only be used with model-fitting functions that require a \code{formula} argument.
#'
#' In contrast to imputation models with categorical exposures, additional arguments need to be specified if the exposure is continuous.
#' All of these additional arguments are related to the sampling procedure for the exposure.
#'
#' Whereas the number of replications \code{nRep} for categorical variables equals the number of levels for the exposure coded as a factor (i.e. the number of hypothetical exposure values), the number of desired replications needs to be specified explicitly for continuous exposures.
#' Its default is 5.
#'
#' If \code{xFit} is left unspecified, the hypothetical exposure levels are automatically sampled from a linear model for the exposure, conditional on a linear combination of all covariates.
#' If one wishes to use another model for the exposure, this default model specification can be overruled by referring to a fitted model object in the \code{xFit} argument.
#' Misspecification of this sampling model does not induce bias in the estimated coefficients and standard errors of the natural effect model.
#'
#' The \code{xSampling} argument allows to specify how the hypothetical exposure levels should be sampled from the conditional exposure distribution (which is either entered explicitly using the \code{xFit} argument or fitted automatically as described in the previous paragraph).
#' The \code{"random"} option randomly samples \code{nRep} draws from the exposure distribution, whereas the \code{"quantiles"} option (default) samples \code{nRep} quantiles at equal-sized probability intervals. Only the latter hence yields fixed exposure levels given \code{nRep} and \code{xFit}. \cr\cr
#' In order to guarantee that the entire support of the distribution is being sampled (which might be a concern if \code{nRep} is chosen to be small), the default lower and upper sampled quantiles are the 5th and 95th percentiles.
#' The intermittent quantiles correspond to equal-sized probability intervals. So, for instance, if \code{nRep = 4}, then the sampled quantiles will correspond to probabilities 0.05, 0.35, 0.65 and 0.95.
#' These default 'outer' quantiles can be changed by specifying the \code{percLim} argument accordingly. By specifying \code{percLim = NULL}, the standard quantiles will be sampled (e.g., 0.2, 0.4, 0.6 and 0.8 if \code{nRep = 4}).
#'
#' @seealso \code{\link{neWeight.default}}, \code{\link{expData}}
#' @examples
#' data(UPBdata)
#' 
#' ## example using glm
#' weightData <- neWeight(negaffect ~ att + gender + educ + age, 
#'                        data = UPBdata, nRep = 2)
#' weights <- attr(weightData, "weights")
#' head(cbind(weightData, weights))
#' 
#' ## example using vglm (yielding identical results as with glm)
#' library(VGAM)
#' weightData2 <- neWeight(negaffect ~ att + gender + educ + age, 
#'                         family = gaussianff, data = UPBdata, nRep = 2, FUN = vglm)
#' weights2 <- attr(weightData2, "weights")
#' head(cbind(weightData2, weights2))
#' @export
neWeight.formula <- function (object, family, data, FUN = glm, nRep = 5, xSampling = c("quantiles", 
    "random"), xFit, percLim = c(0.05, 0.95), ...) 
{
    args <- as.list(match.call())[-1L]
    formula <- object
    nMed <- 1
    joint <- TRUE
    if (missing(family)) 
        family <- formals(FUN)$family
    args$object <- do.call(FUN, eval(list(formula = formula, 
        family = family, data = data, ...)))
    call <- substitute(list(formula = formula, family = family, 
        data = data, ...))
    call[[1]] <- substitute(FUN)
    if (isS4(args$object)) 
        args$object@call <- call
    else args$object$call <- call
    args$formula <- args$family <- args$data <- NULL
    expData <- do.call("neWeight.default", args)
    return(expData)
}
