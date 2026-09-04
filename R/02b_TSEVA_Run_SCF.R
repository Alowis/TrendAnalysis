##########################################################################################
############   THIS SCRIPT IS FOR RUNNING 1 SQUARE OF THE DOMAIN ON THE HPC  ############
##########################################################################################

rm(list = ls())

suppressWarnings(suppressMessages(library(ncdf4)))
suppressWarnings(suppressMessages(library(sf)))
suppressWarnings(suppressMessages(library(rnaturalearth)))
suppressWarnings(suppressMessages(library(rnaturalearthdata)))
suppressWarnings(suppressMessages(library(rgeos)))
suppressWarnings(suppressMessages(library(dplyr)))
suppressWarnings(suppressMessages(library(fs)))
suppressWarnings(suppressMessages(library(tsibble)))
suppressWarnings(suppressMessages(library(ggplot2)))
suppressWarnings(suppressMessages(library(scales)))
suppressWarnings(suppressMessages(library(pracma)))
suppressWarnings(suppressMessages(library(lubridate)))
suppressWarnings(suppressMessages(library(xts)))
suppressWarnings(suppressMessages(library(evd)))
suppressWarnings(suppressMessages(library(POT)))
suppressWarnings(suppressMessages(library(RtsEva)))


###########################  FUNCTIONS   ##################################################

tsEvaTransformSeriesToStationaryMMXTrend <- function(timeStamps, series, timeWindow) {
  tserie <- data.frame(timeStamps, series)
  monthly_max <- tserie %>%
    mutate(month = floor_date(timeStamps, "month")) %>%
    group_by(month) %>%
    summarise(max_value = max(series))

  serieb <- series
  tm <- na.omit(match(as.Date(as.character(monthly_max$month)), as.Date(timeStamps)))
  serieb[-tm] <- NA
  rs <- tsEvaDetrendTimeSeries(timeStamps, serieb, timeWindow)
  detrendSeries <- series - rs@trendSeries
  detrendSerie1 <- serieb - rs@trendSeries
  nRunMn <- rs@nRunMn
  varianceSeries <- tsEvaNanRunningVariance(detrendSerie1, nRunMn)
  varianceSeries <- tsEvaNanRunningMean(varianceSeries, ceiling(nRunMn / 2))
  stdDevSeries <- varianceSeries^0.5
  avgStdDev <- mean(stdDevSeries)
  S <- 2
  N <- timeWindow * 4
  stdDevError <- avgStdDev * (2 * S^2 / N^3)^(1 / 4)
  statSeries <- detrendSeries / stdDevSeries
  statSer3Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn3mom
  statSer4Mom <- tsEvaNanRunningStatistics(statSeries, nRunMn)$rn4mom
  statSer3Mom <- tsEvaNanRunningMean(statSer3Mom, ceiling(nRunMn))
  statSer4Mom <- tsEvaNanRunningMean(statSer4Mom, ceiling(nRunMn))
  trendError <- mean(stdDevSeries) / N^0.5
  trasfData <- list(
    runningStatsMulteplicity = nRunMn, stationarySeries = statSeries,
    trendSeries = rs@trendSeries, trendSeriesNonSeasonal = NULL,
    trendError = trendError, stdDevSeries = stdDevSeries,
    stdDevSeriesNonSeasonal = NULL, stdDevError = stdDevError * rep(1, length(stdDevSeries)),
    timeStamps = timeStamps, nonStatSeries = series,
    statSer3Mom = statSer3Mom, statSer4Mom = statSer4Mom
  )
  return(trasfData)
}

TsEvaNs <- function(timeAndSeries, timeWindow, transfType = "trendPeaks", minPeakDistanceInDays = 10,
                    seasonalityVar = NA, minEventsPerYear = -1, gevMaxima = "annual",
                    ciPercentile = 90, gevType = "GEV", evdType = c("GEV", "GPD"),
                    tail = "high", epy = -1, lowdt = 7, trans = NULL, TrendTh = NA, shape_bnd = NA) {
  timeStamps <- as.POSIXct(timeAndSeries[, 1])
  dt1 <- min(diff(timeStamps), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  if (tdim == "seconds") dt <- dt / 3600
  series <- timeAndSeries[, 2]

  if (epy == -1) {
    if (tail == "high") epy <- 3
    if (tail == "low") epy <- 2
  }

  # Sanity check: remove outlier events >100x the second largest
  sanitycheck <- computeAnnualMaxima(timeAndSeries)
  anmax <- sanitycheck$annualMax[order(sanitycheck$annualMax, decreasing = T)]
  aloc <- sanitycheck$annualMaxIndx[order(sanitycheck$annualMax, decreasing = T)][1]
  x <- anmax[1] / anmax[2]
  if (x > 100) {
    message(paste0("biggest event ", x, " times bigger than second biggest"))
    message("removing this event from timeserie, reruning first steps")
    series[(aloc[1] - 50):(aloc[1] + 50)] <- mean(series)
  }

  valid_types <- c("trend", "seasonal", "trendCIPercentile", "seasonalCIPercentile", "trendPeaks")
  if (!transfType %in% valid_types) {
    stop("\nnonStationaryEvaJRCApproach: transfType can be in (trend, seasonal, trendCIPercentile, trendPeaks)")
  }
  if (minPeakDistanceInDays == -1) stop("label parameter minPeakDistanceInDays must be set")

  minEventsPerYear <- 2

  if (tail == "low") {
    start_index <- 1
    indices_to_extract <- seq(from = start_index, to = length(series), by = lowdt / dt)
    series <- series[indices_to_extract]
    timeStamps <- timeStamps[indices_to_extract]
    if (is.na(shape_bnd)) shape_bnd <- c(-1, 0)
    minEventsPerYear <- 1
    if (trans == "rev") {
      series <- -1 * series
    } else if (trans == "inv") {
      series <- 1 / series
    } else if (trans == "lninv") {
      series <- -log(series)
    }
  } else {
    if (is.na(shape_bnd)) shape_bnd <- c(-0.5, 1)
  }

  if (transfType == "trend") {
    message("\nevaluating long term variations of extremes")
    trasfData <- tsEvaTransformSeriesToStationaryTrendOnly(timeStamps, series, timeWindow)
    gevMaxima <- "annual"
    potEventsPerYear <- epy
  } else if (transfType == "trendChange") {
    message("\nevaluating long term variations of extremes and change point detection")
    trasfData <- tsEvaTransformSeriesToStationaryTrendAndChangepts(timeStamps, series, timeWindow)
    gevMaxima <- "annual"
    potEventsPerYear <- epy
  } else if (transfType == "seasonal") {
    message("\nevaluating long term an seasonal variations of extremes")
    trasfData <- tsEvaTransformSeriesToStationaryMultiplicativeSeasonality(timeStamps, series, timeWindow, seasonalityVar = seasonalityVar)
    gevMaxima <- "monthly"
    potEventsPerYear <- 12
    minEventsPerYear <- 12
  } else if (transfType == "trendCIPercentile") {
    if (is.na(ciPercentile)) stop("For trendCIPercentile transformation the label parameter cipercentile is mandatory")
    message(paste0("\nevaluating long term variations of extremes using the ", ciPercentile, "th percentile"))
    trasfData <- tsEvaTransformSeriesToStationaryTrendOnly_ciPercentile(timeStamps, series, timeWindow, ciPercentile)
    gevMaxima <- "annual"
    potEventsPerYear <- epy
  } else if (transfType == "trendPeaks") {
    print(TrendTh)
    message(paste0("\nevaluating long term variations of the peaks"))
    if (is.na(TrendTh)) {
      TrendTh <- try(tsEvaFindTrendThreshold(series, timeStamps, timeWindow), T)
      if (inherits(TrendTh, "try-error") || length(TrendTh) == 0 || all(is.na(TrendTh))) {
        TrendTh <- 0.1
      }
      trasfData <- tsEvaTransformSeriesToStationaryPeakTrend(timeStamps, series, timeWindow, TrendTh)
    } else {
      if (TrendTh == "MMX") {
        trasfData <- tsEvaTransformSeriesToStationaryMMXTrend(timeStamps, series, timeWindow)
        print("using MMX trend")
      } else {
        trasfData <- tsEvaTransformSeriesToStationaryPeakTrend(timeStamps, series, timeWindow, TrendTh)
        c <- 0
        while (length(unique(trasfData$trendSeries)) < 2) {
          trasfData <- tsEvaTransformSeriesToStationaryPeakTrend(timeStamps, series, timeWindow, TrendTh = TrendTh - c)
          c <- c + 0.1
        }
      }
    }
    gevMaxima <- "annual"
    potEventsPerYear <- epy
  } else if (transfType == "trendChangeCIPercentile") {
    if (is.na(ciPercentile)) stop("For trendCIPercentile transformation the label parameter cipercentile is mandatory")
    message("\n evaluating long term variations of extremes using the ", ciPercentile, "th percentile and change point detection")
    trasfData <- tsEvaTransformSeriesToStationaryTrendAndChangepts_ciPercentile(timeStamps, series, timeWindow, ciPercentile)
    gevMaxima <- "annual"
    potEventsPerYear <- epy
    minEventsPerYear <- 0
  } else if (transfType == "seasonalCIPercentile") {
    if (is.na(ciPercentile)) stop("For seasonalCIPercentile transformation the label parameter cipercentile is mandatory")
    message(paste0("\nevaluating long term variations of extremes using the ", ciPercentile, "th percentile\n"))
    trasfData <- tsEvaTransformSeriesToStatSeasonal_ciPercentile(timeStamps, series, timeWindow, ciPercentile)
    gevMaxima <- "monthly"
    potEventsPerYear <- 12
    minEventsPerYear <- 6
  }

  dtn <- min(diff(trasfData$timeStamps), na.rm = T)
  dtn <- as.numeric(dtn)
  if (dtn < 1) {
    pace <- 1 / dtn
    tsDaily <- seq(1, length(trasfData$timeStamps), by = pace)
    trasfData$stdDevSeriesOr <- trasfData$stdDevSeries
    trasfData$trendSeriesOr <- trasfData$trendSeries
    trasfData$stdDevErrorOr <- trasfData$stdDevError
    trasfData$stdDevSeries <- trasfData$stdDevSeries[tsDaily]
    trasfData$trendSeries <- trasfData$trendSeries[tsDaily]
    trasfData$stdDevError <- trasfData$stdDevError[tsDaily]
  }

  ms <- data.frame(trasfData$timeStamps, trasfData$stationarySeries)
  minPeakDistance <- minPeakDistanceInDays / dtn

  message("\nExecuting stationary eva")
  pointData <- tsEvaSampleDataX(ms, potEventsPerYear, minEventsPerYear, minPeakDistanceInDays, tail, shape_bnd)
  evaAlphaCI <- 0.68
  eva <- tsEVstatisticsX(pointData, evaAlphaCI, gevMaxima, gevType, evdType, shape_bnd)

  if (eva$isValid == FALSE) message("problem in the computation of EVA statistics")

  eva[[2]]$GPDstat$thresholdError <- pointData$POT$thresholdError

  if (eva[[2]][[1]]$method[1] != "No fit") {
    epsilonGevX <- eva[[2]][[1]]$parameters[3]
    errEpsilonX <- epsilonGevX - eva[[2]][[1]]$paramCIs[1, 3]
    sigmaGevX <- eva[[2]][[1]]$parameters[2]
    errSigmaGevX <- sigmaGevX - eva[[2]][[1]]$paramCIs[1, 2]
    muGevX <- eva[[2]][[1]]$parameters[1]
    errMuGevX <- muGevX - eva[[2]][[1]]$paramCIs[1, 1]

    message("\nTransforming to non stationary eva ...\n")
    epsilonGevNS <- epsilonGevX
    errEpsilonGevNS <- errEpsilonX
    sigmaGevNS <- trasfData$stdDevSeries * sigmaGevX
    errSigmaGevFit <- trasfData$stdDevSeries * errSigmaGevX
    errSigmaGevTransf <- sigmaGevX * trasfData$stdDevError
    errSigmaGevNS <- (errSigmaGevTransf^2 + errSigmaGevFit^2)^0.5
    muGevNS <- trasfData$stdDevSeries * muGevX + trasfData$trendSeries
    errMuGevFit <- trasfData$stdDevSeries * errMuGevX
    errMuGevTransf <- ((muGevX * trasfData$stdDevError)^2 + trasfData$trendError^2)^0.5
    errMuGevNS <- (errMuGevTransf^2 + errMuGevFit^2)^0.5

    gevParams <- list(
      epsilon = epsilonGevNS, sigma = sigmaGevNS, mu = muGevNS,
      annualMax = trasfData$nonStatSeries[pointData$annualMaxIndx],
      monthlyMax = trasfData$nonStatSeries[pointData$monthlyMaxIndx],
      annualMaxIndx = pointData$annualMaxIndx,
      monthlyMaxIndx = pointData$monthlyMaxIndx
    )
    if (tolower(gevMaxima) == "annual") {
      gevParams$timeDelta <- 365.25
      gevParams$timeDeltaYears <- 1
    } else if (tolower(gevMaxima) == "monthly") {
      gevParams$timeDelta <- 365.25 / 12
      gevParams$timeDeltaYears <- 1 / 12
    }
    gevParamStdErr <- list(
      epsilonErr = errEpsilonGevNS,
      sigmaErrFit = errSigmaGevFit, sigmaErrTransf = errSigmaGevTransf, sigmaErr = errSigmaGevNS,
      muErrFit = errMuGevFit, muErrTransf = errMuGevTransf, muErr = errMuGevNS
    )
    gevObj <- list(
      method = eva[[2]][[1]]$method, parameters = gevParams,
      paramErr = gevParamStdErr, stationaryParams = eva[[2]][[1]],
      objs = list(monthlyMaxIndexes = pointData$monthlyMaxIndexes)
    )
  } else {
    epsilonGevNS <- eva[[2]][[1]]$parameters[3]
    sigmaGevNS <- trasfData$stdDevSeries * eva[[2]][[1]]$parameters[2]
    muGevNS <- trasfData$stdDevSeries * eva[[2]][[1]]$parameters[1] + trasfData$trendSeries

    gevParams <- list(
      epsilon = epsilonGevNS, sigma = sigmaGevNS, mu = muGevNS,
      annualMax = trasfData$nonStatSeries[pointData$annualMaxIndx],
      monthlyMax = trasfData$nonStatSeries[pointData$monthlyMaxIndx],
      annualMaxIndx = pointData$annualMaxIndx,
      monthlyMaxIndx = pointData$monthlyMaxIndx
    )
    if (tolower(gevMaxima) == "annual") {
      gevParams$timeDelta <- 365.25
      gevParams$timeDeltaYears <- 1
    } else if (tolower(gevMaxima) == "monthly") {
      gevParams$timeDelta <- 365.25 / 12
      gevParams$timeDeltaYears <- 1 / 12
    }
    gevObj <- list(method = "No fit", parameters = gevParams, paramErr = NULL, stationaryParams = NULL, objs.monthlyMaxIndexes = NULL)
  }

  if (eva[[2]][[2]]$method != "No fit") {
    epsilonPotX <- eva[[2]][[2]]$parameters[2]
    errEpsilonPotX <- epsilonPotX - eva[[2]][[2]]$paramCIs[1, 1]
    sigmaPotX <- eva[[2]][[2]]$parameters[1]
    errSigmaPotX <- sigmaPotX - eva[[2]][[2]]$paramCIs[1, 2]
    thresholdPotX <- eva[[2]][[2]]$parameters[3]
    errThresholdPotX <- eva[[2]][[2]]$thresholdError
    nPotPeaks <- eva[[2]][[2]]$parameters[5]
    percentilePotX <- eva[[2]][[2]]$parameters[6]
    dtPeaks <- minPeakDistance
    timeStamps <- as.Date(timeStamps)
    dtPotX <- as.numeric(timeStamps[length(timeStamps)] - timeStamps[1]) / length(series) * dtPeaks

    epsilonPotNS <- epsilonPotX
    errEpsilonPotNS <- errEpsilonPotX
    sigmaPotNS <- sigmaPotX * trasfData$stdDevSeries
    errSigmaPotFit <- trasfData$stdDevSeries * errSigmaPotX
    errSigmaPotTransf <- sigmaPotX * trasfData$stdDevError
    errSigmaPotNS <- (errSigmaPotTransf^2 + errSigmaPotFit^2)^0.5
    thresholdPotNS <- thresholdPotX * trasfData$stdDevSeries + trasfData$trendSeries
    thresholdErrFit <- 0
    thresholdErrTransf <- ((trasfData$stdDevSeries * errThresholdPotX)^2 + (thresholdPotX * trasfData$stdDevError)^2 + trasfData$trendError^2)^0.5
    thresholdErr <- thresholdErrTransf

    potParams <- list(
      epsilon = epsilonPotNS, sigma = sigmaPotNS, threshold = thresholdPotNS,
      percentile = percentilePotX, timeDelta = dtPotX, timeDeltaYears = dtPotX / 365.25,
      timeHorizonStart = min(trasfData$timeStamps), timeHorizonEnd = max(trasfData$timeStamps),
      peaks = trasfData$nonStatSeries[pointData$POT$ipeaks],
      peakID = pointData$POT$ipeaks, peakST = pointData$POT$stpeaks, peakEN = pointData$POT$endpeaks,
      nPeaks = nPotPeaks
    )
    potParamStdErr <- list(
      epsilonErr = errEpsilonPotNS,
      sigmaErrFit = errSigmaPotFit, sigmaErrTransf = errSigmaPotTransf, sigmaErr = errSigmaPotNS,
      thresholdErrFit = thresholdErrFit, thresholdErrTransf = thresholdErrTransf, thresholdErr = thresholdErr
    )
    potObj <- list(method = eva[[2]][[2]]$method, parameters = potParams, paramErr = potParamStdErr, stationaryParams = eva[[2]][[2]], objs = NULL)
  } else {
    dtPeaks <- minPeakDistance
    timeStamps <- as.Date(timeStamps)
    dtPotX <- as.numeric(timeStamps[length(timeStamps)] - timeStamps[1]) / length(series) * dtPeaks
    thresholdPotX <- pointData$POT$threshold
    epsilonPotNS <- pointData$POT$pars[2]
    sigmaPotNS <- pointData$POT$pars[1] * trasfData$stdDevSeries
    thresholdPotNS <- thresholdPotX * trasfData$stdDevSeries + trasfData$trendSeries

    potParams <- list(
      epsilon = epsilonPotNS, sigma = sigmaPotNS, threshold = thresholdPotNS,
      percentile = pointData$POT$percentile, timeDelta = dtPotX, timeDeltaYears = dtPotX / 365.2425,
      timeHorizonStart = min(trasfData$timeStamps), timeHorizonEnd = max(trasfData$timeStamps),
      peaks = trasfData$nonStatSeries[pointData$POT$ipeaks],
      peakID = pointData$POT$ipeaks, peakST = pointData$POT$stpeaks, peakEN = pointData$POT$endpeaks,
      nPeaks = length(pointData$POT$peaks)
    )
    potObj <- list(method = "No fit", parameters = potParams, paramErr = NULL, stationaryParams = NULL, objs = NULL)
  }

  nonStationaryEvaParams <- list(gevObj = gevObj, potObj = potObj)
  return(list(nonStationaryEvaParams = nonStationaryEvaParams, stationaryTransformData = trasfData))
}

tsEVstatisticsX <- function(pointData, alphaCI = 0.95, gevMaxima = "annual", gevType = "GEV", evdType = c("GEV", "GPD"), shape_bnd = c(-0.5, 1)) {
  EVmeta <- list()
  EVdata <- list()
  isValid <- TRUE
  minGEVSample <- 7
  if (is.null(alphaCI)) alphaCI <- 0.95
  if (is.null(gevMaxima)) gevMaxima <- "annual"
  if (is.null(gevType)) gevType <- "GEV"
  if (is.null(evdType)) evdType <- c("GEV", "GPD")

  Tr <- c(5, 10, 20, 50, 100, 200, 500, 1000)
  EVmeta$Tr <- Tr

  methodname <- "GEVstat"
  paramEsts <- numeric(3)
  paramCIs <- matrix(NA, nrow = 2, ncol = 3)
  rlvls <- numeric(length(Tr))

  if (("GEV" %in% evdType) && !is.null(pointData$annualMax)) {
    tmpmat <- if (gevMaxima == "annual") pointData$annualMax else if (gevMaxima == "monthly") pointData$monthlyMax else stop(paste0("Invalid gevMaxima type: ", gevMaxima))
    if (length(tmpmat) >= minGEVSample) {
      tmp <- data.frame(yr = year(pointData$annualMaxDate), dt = tmpmat)
      if (gevType == "GEV") {
        stdfit <- TRUE
        fit <- suppressWarnings(try(evd::fgev(x = tmp$dt, method = "L-BFGS-B", lower = c(-Inf, -Inf, shape_bnd[1]), upper = c(Inf, Inf, shape_bnd[2]), std.err = T), TRUE))
        if (inherits(fit, "try-error")) {
          stdfit <- FALSE
          message("Not able to fit GEV stderr")
          fit <- suppressWarnings(try(evd::fgev(x = tmp$dt, method = "L-BFGS-B", lower = c(-Inf, -Inf, shape_bnd[1]), upper = c(Inf, Inf, shape_bnd[2]), std.err = F), TRUE))
        }
        if (inherits(fit, "try-error")) {
          stdfit <- FALSE
          message("Not able to fit GEV with constrained parameters")
          print(tmp)
          fit <- suppressWarnings(try(evd::fgev(x = tmp$dt, std.err = F), TRUE))
          print(fit)
        }
        paramEsts <- c(mu = fit$par[1], sigma = fit$par[2], xi = fit$par[3])
        alphaCIx <- 1 - alphaCI
        if (stdfit) {
          probs <- c(alphaCIx / 2, 1 - alphaCIx / 2)
          kci <- try(stats::qnorm(probs, paramEsts[3], fit$std.err[3]), TRUE)
          kci[kci < -1] <- -1
          lnsigci <- try(stats::qnorm(probs, log(paramEsts[2]), fit$std.err[2] / paramEsts[2]), silent = T)
          muci <- stats::qnorm(probs, paramEsts[1], fit$std.err[1])
          paramCIs <- cbind(muci = muci, sigci = exp(lnsigci), kci = kci)
        } else {
          paramCIs <- cbind(muci = NA, sigci = NA, kci = NA)
        }
      } else if (gevType == "gumbel") {
        fit <- texmex::evm(y = .data$dt, data = tmp, family = gumbel)
        paramEsts <- c(fit$par[1], exp(fit$par[2]), 0)
        alphaCIx <- 1 - alphaCI
        probs <- c(alphaCIx / 2, 1 - alphaCIx / 2)
        if (is.character(fit$se[1])) {
          message("Gumbel fitted")
        } else {
          message("Gumbel fitted")
          kci <- c(NA, NA)
          lnsigci <- try(stats::qnorm(probs, log(paramEsts[2]), fit$se[2]))
          muci <- stats::qnorm(probs, paramEsts[1], fit$se[1])
          paramCIs <- cbind(muci = muci, sigci = exp(lnsigci), kci = kci)
        }
      } else {
        stop("tsEVstatistics: invalid gevType: ", gevType, ". Can be only GEV, or Gumbel")
      }
      rlvls <- qgev(1 - 1 / Tr, paramEsts[1], paramEsts[2], paramEsts[3])
    } else {
      message("Could not fit GEV")
      methodname <- "No fit"
      rlvls <- NA
      paramEsts <- NA
      paramCIs <- NA
    }
  }
  EVdata$GEVstat <- list(method = methodname, values = rlvls, parameters = paramEsts, paramCIs = paramCIs)

  methodname <- "GPDstat"
  paramEsts <- numeric(6)
  paramCIs <- matrix(NA, nrow = 2, ncol = 3)
  rlvls <- numeric(length(Tr))

  if (("GPD" %in% evdType) && !is.null(pointData$annualMax)) {
    message("Fitted GPD")
    th <- pointData$POT$threshold
    d1 <- pointData$POT$peaks
    fit <- pointData$POT$pars
    if (!inherits(fit, "try-error")) {
      ksi <- fit$param[2]
      sgm <- fit$param[1]
      alphaCIx <- 1 - alphaCI
      probs <- c(alphaCIx / 2, 1 - alphaCIx / 2)
      kci <- try(stats::qnorm(probs, ksi, fit$std.err[2]), silent = T)
      kci[kci < -1] <- -1
      lnsigci <- try(stats::qnorm(probs, log(sgm), fit$std.err[1] / sgm))
      paramCIs <- cbind(kci, sigci = exp(lnsigci))
      paramEstsall <- c(sgm, ksi, pointData$POT$threshold, length(d1), length(pointData$POT$peaks), pointData$POT$percentile)
      rlvls <- pointData$POT$threshold + (sgm / ksi) * ((((length(d1) / length(pointData$POT$peaks)) * (1 / Tr))^(-ksi)) - 1)
      EVdata$GPDstat <- list(method = methodname, values = rlvls, parameters = paramEstsall, paramCIs = paramCIs)
    } else {
      methodname <- "No fit"
      paramEstsall <- c(pointData$POT$pars[1], pointData$POT$pars[2], pointData$POT$threshold, length(d1), length(pointData$POT$peaks), pointData$POT$percentile)
      EVdata$GPDstat <- list(method = methodname, values = NA, parameters = paramEstsall, paramCIs = NA)
      message("could not estimate GPD: bounded distribution")
    }
  } else {
    methodname <- "No fit"
    d1 <- pointData$POT$peaks
    paramEstsall <- c(pointData$POT$pars[1], pointData$POT$pars[2], pointData$POT$threshold, length(d1), length(pointData$POT$peaks), pointData$POT$percentile)
    EVdata$GPDstat <- list(method = methodname, values = NA, parameters = paramEstsall, paramCIs = NA)
    message("could not estimate GPD: bounded distribution")
  }

  return(list(EVmeta = EVmeta, EVdata = EVdata, isValid = isValid))
}

tsEvaSampleDataX <- function(ms, meanEventsPerYear, minEventsPerYear, minPeakDistanceInDays, tail = NA, shape_bnd) {
  pctsDesired <- c(90, 95, 99, 99.9)
  potPercentiles <- c(seq(70, 90, by = 1), seq(91, 95, by = 0.5), seq(95.1, 99.5, by = 0.1))
  if (is.na(tail)) stop("tail for POT selection needs to be 'high' or 'low'")

  POTData <- tsGetPOTX(ms, potPercentiles, meanEventsPerYear, minEventsPerYear, minPeakDistanceInDays, tail, shape_bnd)
  vals <- quantile(ms[, 2], pctsDesired / 100, na.rm = T)
  percentiles <- list(precentiles = pctsDesired, values = vals)

  pointData <- list()
  pointData$completeSeries <- ms
  pointData$POT <- POTData
  pointDataA <- computeAnnualMaxima(ms)
  pointDataM <- computeMonthlyMaxima(ms)
  yrs <- unique(as.numeric(format(as.Date(ms[, 1] + 3600), "%Y")))
  yrs <- yrs - min(yrs)
  pointData$years <- seq(min(yrs), max(yrs), 1)
  pointData$Percentiles <- percentiles
  pointData$annualMax <- pointDataA$annualMax
  pointData$annualMaxDate <- pointDataA$annualMaxDate
  pointData$annualMaxIndx <- pointDataA$annualMaxIndx
  pointData$monthlyMax <- pointDataM$monthlyMax
  pointData$monthlyMaxDate <- pointDataM$monthlyMaxDate
  pointData$monthlyMaxIndx <- pointDataM$monthlyMaxIndx
  return(pointData)
}

tsGetPOTX <- function(ms, pcts, desiredEventsPerYear, minEventsPerYear, minPeakDistanceInDays, tail, shape_bnd) {
  if (minPeakDistanceInDays == -1) stop("label parameter 'minPeakDistanceInDays' must be set")
  dt1 <- min(diff(ms[, 1]), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  if (tdim == "seconds") dt <- dt / 3600
  minPeakDistance <- minPeakDistanceInDays / dt
  minRunDistance <- minPeakDistance
  nyears <- round(as.numeric((max(ms[, 1]) - min(ms[, 1])) / 365.25))
  if (length(pcts) == 1) {
    pcts <- c(pcts - 3, pcts)
    desiredEventsPerYear <- -1
  }

  n_pcts <- length(pcts)
  numperyear <- rep(NA, n_pcts)
  thrsdts <- rep(NA, n_pcts)
  gpp <- rep(NA, n_pcts)
  devpp <- rep(NA, n_pcts)
  fitlist <- vector(mode = "list", length = n_pcts)

  ms_clean <- ms[, 2]
  ms_clean2 <- ms_clean
  if (length(which(is.na(ms_clean))) > 0) {
    ms_clean2 <- ms_clean[-which(is.na(ms_clean))]
    ms_clean[which(is.na(ms_clean))] <- -9999
  }
  if (length(which(is.infinite(ms_clean))) > 0) {
    print("bg")
    ms_clean2 <- ms_clean[-which(is.infinite(ms_clean))]
    ms_clean[which(is.infinite(ms_clean))] <- -9999
  }

  skip <- 0
  for (ipp in 1:n_pcts) {
    if (skip > 0) {
      skip <- skip - 1
      next
    }
    thrsdt <- stats::quantile(ms_clean2, pcts[ipp] / 100, na.rm = TRUE)
    thrsdts[ipp] <- thrsdt
    pks <- if (tail == "high") {
      pracma::findpeaks(ms_clean, minpeakdistance = minPeakDistance, minpeakheight = thrsdt)
    } else {
      declustpeaks(data = ms_clean, minpeakdistance = minPeakDistance, minrundistance = minRunDistance, qt = thrsdt)
    }
    if (is.null(pks) || length(pks) == 0) next
    numperyear[ipp] <- nrow(pks) / nyears
    penalty_strength <- 1000
    penalty_factor <- 1
    if (numperyear[ipp] < minEventsPerYear) {
      deficit <- 1 - (numperyear[ipp] / minEventsPerYear)
      penalty_factor <- 1 + (deficit^2 * penalty_strength)
    }
    if (numperyear[ipp] <= (desiredEventsPerYear + 1)) {
      fgpd <- suppressWarnings(try(POT::fitgpd(pks[, 1],
        threshold = thrsdt, est = "mle",
        method = "L-BFGS-B", lower = c(1e-6, shape_bnd[1]), upper = c(Inf, shape_bnd[2]),
        std.err.type = "observed"
      ), silent = TRUE))
      fitlist[[ipp]] <- fgpd
      if (inherits(fgpd, "try-error")) {
        devpp[ipp] <- 1e9
        gpp[ipp] <- NA
      } else {
        devpp[ipp] <- fgpd$deviance + penalty_factor
        gpp[ipp] <- fgpd$fitted.values[2]
      }
    }
  }

  pkx <- declustpeaks(data = ms_clean, minpeakdistance = minPeakDistance, minrundistance = minRunDistance, qt = stats::quantile(ms_clean, pcts[1] / 100, na.rm = T))
  devpp[1] <- NA
  dist_to_bnd <- pmin(abs(gpp - shape_bnd[1]), abs(gpp - shape_bnd[2]))
  dshap <- c(0, diff(gpp))
  low_tail_threshold <- shape_bnd[2]
  devpx <- devpp
  if (tail == "low") devpx[gpp >= low_tail_threshold] <- NA
  devpx[abs(dshap) > 0.5] <- NA
  print(devpx)

  trip_in_bounds <- which(!is.na(devpx) & gpp >= shape_bnd[1] & gpp <= shape_bnd[2])
  if (length(trip_in_bounds) > 0) {
    trip <- trip_in_bounds[which.min(devpx[trip_in_bounds])]
    message(paste0("Optimal shape found: ", round(gpp[trip], 3)))
  } else {
    valid_indices <- which(!is.na(devpx))
    print(valid_indices)
    if (length(valid_indices) == 0) {
      valid_indices <- 1:length(devpp)
      devpx <- devpp
    }
    v_dev <- devpx[valid_indices]
    v_dist <- dist_to_bnd[valid_indices]
    norm_dev <- (v_dev - min(v_dev, na.rm = T)) / (diff(range(v_dev, na.rm = T)) + 1e-10)
    norm_dist <- (v_dist - min(v_dist, na.rm = T)) / (diff(range(v_dist, na.rm = T)) + 1e-10)
    combined_score <- (0.5 * norm_dev) + (0.5 * norm_dist)
    trip <- valid_indices[which.min(combined_score)]
    message(paste0("Shape outside boundaries. Fallback selected: ", round(gpp[trip], 3)))
  }

  print(numperyear)
  message(paste0("\nmax threshold is: ", pcts[trip], "%"))
  message(paste0("\nshape parameter is: ", round(gpp[trip], 2)))
  message(paste0("\naverage number of events per year = ", round(numperyear[trip], 1)))

  diffNPerYear <- mean(diff(stats::na.omit(rev(numperyear)), na.rm = TRUE))
  if (diffNPerYear == 0) diffNPerYear <- 1
  diffNPerYear <- 1
  thresholdError <- -mean(diff(stats::na.omit(thrsdts)) / diffNPerYear) / 2
  indexp <- trip
  print(indexp)
  if (length(indexp) > 0) {
    thrsd <- stats::quantile(ms_clean2, pcts[indexp] / 100)
    pct <- pcts[indexp]
    fit <- fitlist[[indexp]]
  } else {
    thrsd <- 0
    pct <- 0
  }

  if (tail == "high") pks_and_locs <- pracma::findpeaks(ms_clean, minpeakdistance = minPeakDistance, minpeakheight = thrsd)
  if (tail == "low") pks_and_locs <- declustpeaks(data = ms_clean, minpeakdistance = minPeakDistance, minrundistance = minRunDistance, qt = thrsd)

  POTdata <- list(
    threshold = thrsd, thresholdError = thresholdError, percentile = pct,
    peaks = pks_and_locs[, 1], stpeaks = pks_and_locs[, 3], endpeaks = pks_and_locs[, 4],
    ipeaks = pks_and_locs[, 2], time = ms[pks_and_locs[, 2], 1], pars = fit
  )
  return(POTdata)
}

outletopen <- function(dir, outletname, nrspace = rep(NA, 5)) {
  ncbassin <- paste0(dir, "/", outletname, ".nc")
  ncb <- nc_open(ncbassin)
  namev <- names(ncb[["var"]])[1]
  if (!is.na(nrspace[1])) {
    start <- as.numeric(nrspace[c(2, 4)])
    count <- as.numeric(nrspace[c(3, 5)]) - start + 1
  } else {
    start <- c(1, 1)
    count <- c(llo, lla)
  }
  londat <- ncvar_get(ncb, "lon", start = start[1], count = count[1])
  llo <- length(londat)
  latdat <- ncvar_get(ncb, "lat", start = start[2], count = count[2])
  lla <- length(latdat)
  outlets <- as.vector(ncvar_get(ncb, namev, start = start, count = count))
  outll <- expand.grid(londat, latdat)
  lonlatloop <- expand.grid(c(1:llo), c(1:lla))
  outll$idlo <- lonlatloop$Var1
  outll$idla <- lonlatloop$Var2
  outll <- outll[which(!is.na(outlets)), ]
  outlets <- outlets[which(!is.na(outlets))]
  return(data.frame(outlets, outll))
}

disNcopenloc <- function(fname, dir, outloc, idc) {
  ncd <- nc_open(paste0(dir, "/", fname, ".nc"))
  namev <- names(ncd[["var"]])[1]
  time <- ncvar_get(ncd, "time")
  lt <- length(time)
  londat <- ncvar_get(ncd, "lon")
  latdat <- ncvar_get(ncd, "lat")
  start <- c(outloc$idlo[idc], outloc$idla[idc], 1)
  count <- c(1, 1, lt)
  outlets <- as.vector(ncvar_get(ncd, namev, start = start, count = count))
  outll <- data.frame(
    outlets,
    outid = rep(outloc[idc, 1], lt),
    lon = rep(londat[idc], lt),
    lat = rep(latdat[idc], lt),
    time
  )
  return(outll)
}

ComputeReturnLevels <- function(nonStationaryEvaParams, RPgoal, timeIndex) {
  epsilonGEV <- nonStationaryEvaParams[[1]]$parameters$epsilon
  sigmaGEV <- mean(nonStationaryEvaParams[[1]]$parameters$sigma[timeIndex])
  muGEV <- mean(nonStationaryEvaParams[[1]]$parameters$mu[timeIndex])
  epsilonGPD <- nonStationaryEvaParams[[2]]$parameters$epsilon
  sigmaGPD <- mean(nonStationaryEvaParams[[2]]$parameters$sigma[timeIndex])
  thresholdGPD <- mean(nonStationaryEvaParams[[2]]$parameters$threshold[timeIndex])
  nPeaks <- nonStationaryEvaParams[[2]]$parameters$nPeaks
  thStart <- nonStationaryEvaParams[[2]]$parameters$timeHorizonStart
  thEnd <- nonStationaryEvaParams[[2]]$parameters$timeHorizonEnd
  sampleTimeHorizon <- as.numeric((thEnd - thStart) / 365.2425)

  if (nonStationaryEvaParams[[1]]$method == "No fit") {
    print("could not fit EVD to this pixel")
    ParamGEV <- c(epsilonGEV, sigmaGEV, muGEV, NA, NA, NA)
    names(ParamGEV) <- c("epsilonGEV", "sigmaGEV", "muGEV", "epsilonStdErrGEV", "sigmaStdErrGEV", "muStdErrGEV")
    ParamGPD <- c(epsilonGPD, sigmaGPD, thresholdGPD, NA, NA, NA, nPeaks, sampleTimeHorizon)
    names(ParamGPD) <- c("epsilonGPD", "sigmaGPD", "thresholdGPD", "epsilonStdErrGPD", "sigmaStdErrGPD", "thresholdStdErrGPD", "nPeaks", "SampleTimeHorizon")
    return(list(Fit = "No fit", Params = c(ParamGEV, ParamGPD)))
  }

  epsilonStdErrGEV <- nonStationaryEvaParams[[1]]$paramErr$epsilonErr
  sigmaStdErrGEV <- mean(nonStationaryEvaParams[[1]]$paramErr$sigmaErr[timeIndex])
  muStdErrGEV <- mean(nonStationaryEvaParams[[1]]$paramErr$muErr[timeIndex])
  epsilonStdErrGPD <- nonStationaryEvaParams[[2]]$paramErr$epsilonErr
  sigmaStdErrGPD <- mean(nonStationaryEvaParams[[2]]$paramErr$sigmaErr[timeIndex])
  thresholdStdErrGPD <- mean(nonStationaryEvaParams[[2]]$paramErr$thresholdErr[timeIndex])

  returnLevelsGEV <- tsEvaComputeReturnLevelsGEV(epsilonGEV, sigmaGEV, muGEV, epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV, RPgoal)
  returnLevelsGPD <- tsEvaComputeReturnLevelsGPD(epsilonGPD, sigmaGPD, thresholdGPD, epsilonStdErrGPD, sigmaStdErrGPD, thresholdStdErrGPD, nPeaks, sampleTimeHorizon, RPgoal)

  ParamGEV <- c(epsilonGEV, sigmaGEV, muGEV, epsilonStdErrGEV, sigmaStdErrGEV, muStdErrGEV)
  names(ParamGEV) <- c("epsilonGEV", "sigmaGEV", "muGEV", "epsilonStdErrGEV", "sigmaStdErrGEV", "muStdErrGEV")
  ParamGPD <- c(epsilonGPD, sigmaGPD, thresholdGPD, epsilonStdErrGPD, sigmaStdErrGPD, thresholdStdErrGPD, nPeaks, sampleTimeHorizon)
  names(ParamGPD) <- c("epsilonGPD", "sigmaGPD", "thresholdGPD", "epsilonStdErrGPD", "sigmaStdErrGPD", "thresholdStdErrGPD", "nPeaks", "SampleTimeHorizon")

  return(list(
    Fit = "Fitted",
    ReturnLevels = c(
      ReturnPeriod = RPgoal,
      GEV = as.numeric(returnLevelsGEV$returnLevels), GPD = as.numeric(returnLevelsGPD$returnLevels),
      errGEV = as.numeric(returnLevelsGEV$returnLevelsErr), errGPD = as.numeric(returnLevelsGPD$returnLevelsErr)
    ),
    Params = c(ParamGEV, ParamGPD)
  ))
}

RPcalc <- function(params, RPiGEV, RPiGPD) {
  paramx <- data.frame(t(params))
  qxV <- 1 - exp(-(1 + paramx$epsilonGEV * (RPiGEV - paramx$muGEV) / paramx$sigmaGEV)^(-1 / paramx$epsilonGEV))
  returnPeriodGEV <- if (is.na(qxV)) 9999 else 1 / qxV
  X0 <- paramx$nPeaks / paramx$SampleTimeHorizon
  qxD <- (1 + paramx$epsilonGPD * (RPiGPD - paramx$thresholdGPD) / paramx$sigmaGPD)^(-1 / paramx$epsilonGPD)
  returnPeriodGPD <- if (is.na(qxD)) 9999 else 1 / (X0 * qxD)
  return(c(GEV = returnPeriodGEV, GPD = returnPeriodGPD))
}

interid <- function(data, trans, WindowSize) {
  dt1 <- min(diff(data$date), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  nRunMn <- ceiling(WindowSize / dt)
  data$Q7 <- rollmean(data$Qs, nRunMn, align = "right", fill = NA)
  if (trans == "rev") {
    data$Qtrans <- -(data$Q7)
  } else if (trans == "inv") {
    data$Qtrans <- 1 / data$Q7
  } else if (trans == "lninv") {
    data$Qtrans <- -ln(data$Q7)
  }
  if (length(which(is.na(data$Q7))) == length(data$Q7)) {
    print("no data in this pixel")
    list0 <- NA
    dis07 <- data
    l0 <- NA
    fl <- NA
    mindis <- NA
  } else {
    if (length(which(is.na(data$Qs))) > 0) print("Na alert")
    mindis <- min(data$Q7, na.rm = T)
    m0 <- length(which(data$Qs7 == mindis))
    start_index <- 1
    indices_to_extract <- seq(from = start_index, to = length(data$Q7), by = WindowSize / dt)
    datat <- data$Q7[indices_to_extract]
    l0 <- length(which(datat <= 1e-4))
    dayysbelow <- tsEvaNanRunnigBlowTh(series = data$Q7, threshold = 1e-4, windowSize = 4 * 365 * 30)
    dayysbelow$time <- as.Date(data$date[dayysbelow$time])
    yrtot <- length(unique(year(data$date)))
    list0 <- NA
    fl <- 0
    if (l0 >= 7) {
      print("intermittent river class2")
      fl <- 3
      data$Qinv <- NA
      list0 <- data$date[which(data$Qs == 0)]
    } else if (l0 > 1) {
      print("intermittent river class1 ")
      fl <- 1
      list0 <- data$date[which(data$Q7 == mindis)]
    } else if (mindis > 0 & m0 >= yrtot) {
      print(paste0("river with floor low flow ", mindis))
      fl <- 2
      list0 <- data$date[which(data$Q7 == mindis)]
    }
    dis07 <- data[, c(2, 3, 4)]
  }
  return(list(zerodate = list0, trdis = dis07, DaysBlow = dayysbelow, flags = c(n0d = l0, intertype = fl, mindischarge = mindis)))
}

UpAopen <- function(dir, outletname, Sloc_final) {
  ncb <- nc_open(paste0(dir, outletname))
  namev <- names(ncb[["var"]])[2]
  londat <- ncvar_get(ncb, "lon")
  latdat <- ncvar_get(ncb, "lat")
  llo <- length(londat)
  lla <- length(latdat)
  outlets <- as.vector(ncvar_get(ncb, namev, start = c(1, 1), count = c(llo, lla))) / 1000000
  outll <- expand.grid(londat, latdat)
  lonlatloop <- expand.grid(c(1:llo), c(1:lla))
  outll$upa <- outlets
  outll$idlo <- lonlatloop$Var1
  outll$idla <- lonlatloop$Var2
  outll$latlong <- paste(round(outll$Var1, 4), round(outll$Var2, 4), sep = " ")
  return(inner_join(outll, Sloc_final, by = "latlong"))
}


###########################  ARGUMENTS   #################################################

# --- Path config -------------------------------------------------------------
# This script defines its own functions inline (self-contained for the HPC),
# but the data root is taken from config_paths.R so it stays in sync with the
# rest of the pipeline. Resolve the config location without rstudioapi so it
# also works under Rscript.
.cfg_dir <- {
  a <- commandArgs(FALSE)
  fa <- grep("^--file=", a, value = TRUE)
  if (length(fa) > 0) dirname(normalizePath(sub("^--file=", "", fa[1]))) else getwd()
}
source(file.path(.cfg_dir, "config_paths.R"))
workDir <- hydroDir # canonical data root (no trailing slash; use file.path)

# --- HPC mode: arguments passed via command line ---
# To run on HPC: Rscript script.R <Nsq> <haz> <sce> <startid> <endid>
args <- commandArgs(TRUE)
if (length(args) >= 1) {
  argus <- as.vector(unlist(strsplit(args, split = " ")))
  Nsq <- as.numeric(argus[1])
  haz <- ifelse(length(argus) >= 2, argus[2], "drought")
  sce <- ifelse(length(argus) >= 3, argus[3], "SCF")
  startid <- ifelse(length(argus) >= 4, as.numeric(argus[4]), 1)
  endid <- ifelse(length(argus) >= 5, as.numeric(argus[5]), NA)
} else {
  # --- Local defaults ---
  haz <- "drought"
  Nsq <- 42
  sce <- "SCF"
  startid <- 1
  endid <- NA
}
tail <- "high"
var <- "dis"
outlets <- "RNetwork"
outletname <- "GeoData/efas_rnet_100km_01min"
season <- "nonfrost"


rspace <- read.csv(file.path(workDir, "subspace_efas.csv"))
rspace <- rspace[, -1]
nrspace <- rspace[Nsq, ]
print(nrspace)
outhybas <- outletopen(workDir, outletname, nrspace)
Idstart <- as.numeric(Nsq) * ID_MULT
if (length(outhybas$outlets) > 0) {
  outhybas$outlets <- seq((Idstart + 1), (Idstart + length(outhybas$outlets)))
}
unikout <- outhybas$outlets
outhybas$latlong <- paste(round(outhybas$Var1, 4), round(outhybas$Var2, 4), sep = " ")

UpArea <- UpAopen(workDir, "/GeoData/upArea_European_01min.nc", outhybas)
outhybas$upa <- UpArea$upa
#
# # Load SocCF results
# if (haz == "drought") namefile <- "Drought.nonfrost.SocCF"
# if (haz == "flood") namefile <- "flood.year.SocCF"
# load(file = paste0(workDir, "/", haz, "/peaks.", namefile, ".Rdata"))
# load(file = paste0(workDir, "/", haz, "/params.", namefile, ".Rdata"))
#
# Paramsfl$square <- round(Paramsfl$catchment / 100000)
# Paramsfl <- Paramsfl[which(Paramsfl$square == Nsq), ]
# Peaksave$square <- round(Peaksave$catch / 100000)
# Peaksave <- Peaksave[which(Peaksave$square == Nsq), ]

filename <- file.path("RiverData", paste0("dis_", Nsq, "_1951_2020_scf_RNetwork"))

dists <- disNcopenloc(filename, workDir, outhybas, 1)
df.dis <- dists
print(paste0("hazard: ", haz, " opening square ", Nsq, " /88"))
timeStamps <- as.POSIXct(as.Date(df.dis$time, origin = "1979-01-01") - 1 / 24)
txx <- timeStamps
df.dis$timeStamps <- txx
names(df.dis)[c(1, 2)] <- c("dis", "outlets")

# Load drought-specific data
haz <- "drought"
if (haz == "drought") {
  load(file = file.path(droughtDir, "catchment_frost.Rdata"))
  frostcat <- frostcat[-which(year(frostcat$time) == 1950), ]
  frostcat <- frostcat[-1, ]
  Catchmentrivers7 <- read.csv(file.path(geoDir, "HYBAS07/from_hybas_eu_onlyid.csv"), encoding = "UTF-8", header = T, stringsAsFactors = F)
  outhyb07 <- outletopen(workDir, "GeoData/HYBAS07/outletsv8_hybas07_01min", nrspace)
  mycat <- Catchmentrivers7[match(outhyb07$outlets, Catchmentrivers7$pointid), ]
  hybas07 <- read_sf(dsn = file.path(geoDir, "HYBAS07/hybas_eu_lev07_v1c.shp"))
  Catamere07 <- inner_join(fortify(hybas07), Catchmentrivers7, by = "HYBAS_ID")
  Catamere07$llcoord <- paste(round(Catamere07$POINT_X, 4), round(Catamere07$POINT_Y, 4), sep = " ")
  Catf7 <- inner_join(Catamere07, outhybas, by = c("llcoord" = "latlong"))
  st_geometry(Catf7) <- NULL
  tail <- "low"
}



# Load and merge thresholds (produced by 01_TSEVA_TrendThresholdSel.R).
# 01 writes files named trenTH_x_<sce>_<tail>_<Nsq>.csv using the ID_MULT
# scheme, so no id remapping is needed here.

TH2 <- read.csv(file.path(threshDir, paste0("trenTH_SCF_", tail, "_", Nsq, ".csv")))
TH3 <- TH2

# retain thresholds from SCF run unless it is NA
thresh_vec <- data.frame(TH3$cid, TH3$Th_new)

names(thresh_vec) <- c("cid", "th")
#id conversion for matching
thresh_vec$cid <- as.numeric(thresh_vec$cid - Nsq*10000 + Nsq*100000)
Nsq <- as.numeric(Nsq)

# Loop bounds: honour command-line args, else run the whole square
if (is.na(endid)) endid <- length(unikout)
 endid <- 7   # (debug) limit to first pixels
RetPerGPD <- c()
RetPerGEV <- c()
RetLevGEV <- c()
RetLevGPD <- c()
parlist <- c()
peaklist <- c()
catlist <- c()
IRES <- c()

for (idfix in startid:endid) {
  start_time <- Sys.time()
  print(paste0("hazard:", haz, " square: ", Nsq, " pixel: ", idfix, "/", endid))
  catch <- as.numeric(unikout[idfix])
  upa <- as.numeric(outhybas$upa[idfix])

  timeStamps <- txx
  thresh <- thresh_vec$th[which(thresh_vec$cid == catch)]
  frosttime <- NA

  df.disX <- disNcopenloc(filename, workDir, outhybas, idfix)
  series <- data.frame(txx, df.disX$outlets)
  names(series) <- c("date", "Qs")
  rmv <- which(as.integer(format(series$date, "%Y")) == 1950)
  if (length(rmv) > 0) series <- series[-rmv, ]
  series <- series[-1, ]

  if (haz == "drought") {
    trans <- "rev"
    catmat <- Catf7[which(Catf7$outlets == catch), ]
    Tcatmat <- mycat[which(mycat$HYBAS_ID == catmat$HYBAS_ID), ]
    Tcatchment <- which(colnames(frostcat) == Tcatmat$pointid)
    intermit <- interid(series, trans, WindowSize = 7)
    interflag <- intermit$flags[2]
    series <- data.frame(series$date, intermit$trdis$Q7)
    if (length(Tcatchment) > 0) {
      frostserie <- data.frame(frostcat[, 1], frostcat[, Tcatchment])
      frosttime <- which(frostserie[, 2] < 0)
    } else {
      frosttime <- NA
    }
    ciPercentile <- 80
    minPeakDistanceInDays <- 30
    tail <- "low"
  } else if (haz == "flood") {
    ciPercentile <- 95
    minPeakDistanceInDays <- 5 + log(upa / 2.59) # dynamic peak distance based on upstream area
    interflag <- 0
    series <- max_daily_value(series)
    tail <- "high"
    trans <- "ori"
  }

  names(series) <- c("timestamp", "dis")
  dt1 <- min(diff(series$timestamp), na.rm = T)
  dt <- as.numeric(dt1)
  tdim <- attributes(dt1)$units
  if (tdim == "hours") dt <- dt / 24
  timeDays <- if (dt == 1) series$timestamp else unique(as.Date(series$timestamp))

  bounds <- c(year(timeDays[1]), year(timeDays[length(timeDays)]))
  tbound <- c(as.Date(paste0(bounds[1], "-12-31")), as.Date(paste0(bounds[2], "-12-31")))
  Impdates <- seq(tbound[1], tbound[2], by = "1 year")

  nv <- length(unique(series$dis))
  if (length(na.omit(series$dis)) > 1 & interflag < 3 & nv > 15) {
    if (length(which(is.na(series$dis))) > 0) {
      print("Na alert")
      series$dis <- tsEvaFillSeries(series$timestamp, series$dis)
    }
    timeAndSeries <- series
    names(timeAndSeries) <- c("timestamp", "data")

    if (haz == "drought" & length(frosttime[!is.na(frosttime)]) > 1) {
      if (season == "nonfrost") {
        print("nonfrost season")
        timeAndSeries$data[frosttime] <- NA
      } else if (season == "frost") {
        print("frost season")
        timeAndSeries$data[-frosttime] <- NA
      } else if (season == "year") {
        print("no seasonal divide")
      } else {
        print("season must be frost or nonfrost")
      }
    } else {
      print("no frost season for this river")
    }

    series <- timeAndSeries[, 2]
    timeWindow <- 365.25 * 30
    timeStamps <- timeAndSeries$timestamp
    cat(paste0("\nsquare: ", Nsq, " pixel: ", idfix, "/", endid))

    Nonstat <- TsEvaNs(timeAndSeries, timeWindow,
      transfType = "trendPeaks",
      ciPercentile = ciPercentile, minPeakDistanceInDays = minPeakDistanceInDays,
      lowdt = 7, trans = trans, tail = tail, TrendTh = thresh
    )
    nonStationaryEvaParams <- Nonstat[[1]]
    stationaryTransformData <- Nonstat[[2]]

    stationaryTransformData$timeStampsDay <- unique(as.Date(stationaryTransformData$timeStamps))
    pikos <- data.frame(
      value = nonStationaryEvaParams$potObj$parameters$peaks,
      timeID = nonStationaryEvaParams$potObj$parameters$peakID,
      tIDstart = nonStationaryEvaParams$potObj$parameters$peakST,
      tIDend = nonStationaryEvaParams$potObj$parameters$peakEN
    )
    pikos$time <- timeStamps[pikos$timeID]
    pikos$catch <- catch

    dt1 <- min(diff(timeStamps), na.rm = T)
    dt <- as.numeric(dt1)
    tdim <- attributes(dt1)$units
    if (tdim == "hours") dt <- dt / 24
    timeDays <- if (dt == 1) stationaryTransformData$timeStamps else stationaryTransformData$timeStampsDay

    tbound <- c(as.Date(paste0(year(timeDays[1]), "-12-31")), as.Date(paste0(year(timeDays[length(timeDays)]), "-12-31")))
    Impdates <- seq(tbound[1], tbound[2], by = "1 years")
    dtect <- c(diff(yday(timeDays)), -1)
    last_days <- timeDays[which(dtect < 0)]
    tindexes <- match(last_days, timeDays)

    RPgoal <- 10
    RLevs100 <- ComputeReturnLevels(nonStationaryEvaParams, RPgoal, tindexes[1])

    if (RLevs100$Fit == "No fit") {
      RLgev <- RLgpd <- nonStationaryEvaParams$gevObj$parameters$annualMax
      names(RLgev) <- names(RLgpd) <- year(Impdates)
      nRPgev <- nRPgpd <- rep(NA, length(Impdates))
      names(nRPgev) <- names(nRPgpd) <- year(Impdates)
      params <- data.frame(matrix(ncol = 19, nrow = (length(Impdates) - 1)))
      params[, 1] <- catch
      params[, 2] <- year(Impdates)[-1]
      params[, 5] <- interflag
      colnames(params) <- if (is.null(colnames(parlist))) rep("nom", 19) else colnames(parlist)
    } else {
      RLgev <- RLevs100$ReturnLevels[2]
      RLgpd <- RLevs100$ReturnLevels[3]
      ERgev <- RLevs100$ReturnLevels[4]
      ERgpd <- RLevs100$ReturnLevels[5]
      nRPgev <- nRPgpd <- 10
      params <- c()
      for (t in 2:length(Impdates)) {
        timeIndex <- tindexes[t]
        RLevs100i <- ComputeReturnLevels(nonStationaryEvaParams, RPgoal, timeIndex)
        params <- c(catch, year(Impdates[t]), timeIndex, minPeakDistanceInDays, RLevs100i$Params, nonStationaryEvaParams$potObj$parameters$percentile)
        names(params)[1:4] <- c("catchment", "Year", "timeIndex", "minpeakdistance")
        names(params)[19] <- "percentile"
        Rper <- RPcalc(params, RPiGEV = RLevs100$ReturnLevels[2], RPiGPD = RLevs100$ReturnLevels[3])
        nRPgpd <- c(nRPgpd, Rper[2])
        nRPgev <- c(nRPgev, Rper[1])
        RLgev <- cbind(RLgev, RLevs100i$ReturnLevels[2])
        RLgpd <- cbind(RLgpd, RLevs100i$ReturnLevels[3])
        ERgev <- cbind(ERgev, RLevs100i$ReturnLevels[4])
        ERgpd <- cbind(ERgpd, RLevs100i$ReturnLevels[5])
        if (length(parlist) > 1) colnames(parlist) <- names(params)
        parlist <- rbind(parlist, params)
      }
      RLgev <- as.data.frame(RLgev)
      names(RLgev) <- year(Impdates)
      rownames(RLgev) <- RPgoal
      RLgpd <- as.data.frame(RLgpd)
      names(RLgpd) <- year(Impdates)
      rownames(RLgpd) <- RPgoal
      nRPgev <- as.data.frame(t(nRPgev))
      names(nRPgev) <- year(Impdates)
      nRPgpd <- as.data.frame(t(nRPgpd))
      names(nRPgpd) <- year(Impdates)
      peaklist <- rbind(peaklist, pikos)
    }
  } else {
    cat(paste0("\n No values in this pixel ", idfix, " \n or intermittent river (flag = ", interflag, ")"))
    if (is.na(interflag)) interflag <- -9999
    if (interflag > 0) {
      datex <- yday(intermit$DaysBlow$time)
      dtect <- c(diff(datex), -1)
      last_days <- intermit$DaysBlow$time[which(dtect < 0)]
      tindexes <- match(last_days, intermit$DaysBlow$time)
      oops <- intermit$DaysBlow[tindexes, ]
    } else {
      filling <- NA
    }
    RLgev <- RLgpd <- oops$RP
    names(RLgev) <- names(RLgpd) <- year(Impdates)
    nRPgev <- nRPgpd <- rep(NA, length(Impdates))
    names(nRPgev) <- names(nRPgpd) <- year(Impdates)
    params <- data.frame(matrix(ncol = 19, nrow = length(Impdates) - 1))
    params[, 1] <- catch
    params[, 2] <- year(Impdates)[-1]
    params[, 5] <- interflag
    print(colnames(parlist))
    colnames(params) <- if (is.null(colnames(parlist))) {
      print("hello")
      rep("nom", 19)
    } else {
      colnames(parlist)
    }
    parlist <- as.data.frame(rbind(parlist, params))
    pikos <- data.frame(matrix(ncol = 6, nrow = 1))
    pikos[, 1] <- NA
    pikos[, 6] <- catch
    colnames(pikos) <- if (is.null(colnames(peaklist))) c("value", "timeID", "tIDstart", "tIDend", "time", "catch") else colnames(peaklist)
    print("bbbb")
    peaklist <- as.data.frame(rbind(peaklist, pikos))
  }

  catlist <- c(catlist, catch)
  IRES <- c(IRES, interflag)
  RetLevGEV <- rbind(RetLevGEV, RLgev)
  RetLevGPD <- rbind(RetLevGPD, RLgpd)
  RetPerGEV <- rbind(RetPerGEV, nRPgev)
  RetPerGPD <- rbind(RetPerGPD, nRPgpd)
  cat(paste0("\nloop duration: ", round(Sys.time() - start_time, 2), " seconds\n"))
}

Results <- list(
  parameters = parlist,
  RetLevGEV = RetLevGEV, RetLevGPD = RetLevGPD,
  RetPerGEV = RetPerGEV, RetPerGPD = RetPerGPD,
  Peaks = peaklist, catrest = data.frame(catlist, IRES)
)

hazDir <- file.path(hydroDir, haz)
if (!dir.exists(hazDir)) dir.create(hazDir, recursive = TRUE, showWarnings = FALSE)

save(Results, file = file.path(hazDir, paste0("ResCat6h_", outlets, "_", Nsq, "_1951_2020_revFt.Rdata")))

# for subsequent analysis
save(parlist, file = file.path(hazDir, paste0("Params_", haz, ".", season, "_", Nsq, "_SocCF.Rdata")))
save(parlist, file = file.path(hazDir, paste0("RL100x_", haz, ".", season, "_", Nsq, "_SocCF.Rdata")))
save(parlist, file = file.path(hazDir, paste0("peaks_", haz, ".", season, "_", Nsq, "_SocCF.Rdata")))
