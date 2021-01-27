#' Analyze spatial count data using the flexible spatial scan statistic
#'
#' The rflexscan package provides functions and classes to analyze spatial count
#' data using the flexible spatial scan statistic developed by Tango and
#' Takahashi (2005). This package designed for any of the following interrelated
#'  purposes:
#' \enumerate{
#'   \item To evaluate reported spatial disease clusters, to see if they are
#'         statistically significant.
#'   \item To test whether a disease is randomly distributed over space.
#'   \item To perform geographical surveillance of disease, to detect areas of
#'         significantly high rates.
#' }
#' This package implements a wrapper for the C routine used in the FleXScan 3.1.2
#' developed by Takahashi, Yokoyama, and Tango.
#'
#' @references
#' \itemize{
#'   \item Tango T. and Takahashi K. (2005). A flexibly shaped spatial scan
#'   statistic for detecting clusters, International Journal of Health
#'   Geographics 4:11.
#'   \item Takahashi K, Yokoyama T and Tango T. (2010). FleXScan v3.1: Software
#'   for the Flexible Scan Statistic. National Institute of Public Health, Japan,
#'   \url{https://sites.google.com/site/flexscansoftware/home}.
#' }
#'
#' @aliases NULL multicluster-package
#'
"_PACKAGE"

flexscan.model <- c("POISSON", "BINOMIAL")
flexscan.stattype <- c("ORIGINAL", "RESTRICTED")
flexscan.scanmethod <- c("FLEXIBLE", "CIRCULAR")
flexscan.rantype <- c("MULTINOMIAL", "POISSON")


#' Detect multiple spatial clusters using the flexible/circular scan statistic
#'
#' This function analyzes spatial count data using the flexible spatial scan
#' statistic developed by Tango and Takahashi (2005) or Kulldorff's circular
#' spatial scan statistic (1997), and detect spatial disease clusters.
#'
#' @details
#' Centroid coordinates for each region should be specified EITHER by Cartesian
#' coordinates using arguments \code{x} and \code{y} or by latitudes and
#' longitudes using arguments \code{lat} and \code{lon}.
#' Note that \code{lat} and \code{lon} are DEPRECATED due to accuracy issues.
#' This feature is implemented to maintain compatibility with FleXScan software.
#' We recommend to transform latitude and longitude onto the Cartesian
#' coordinate system beforehand (using \code{spTransform} function in sp package,
#' for example) and use the \code{x} and \code{y} parameters that are projected
#' coordinates.
#'
#' @param x
#' A vector of X-coordinates.
#'
#' @param y
#' A vector of Y-coordinates.
#'
#' @param observed
#' A vector with the observed number of disease cases.
#'
#' @param expected
#' A vector with the expected number of disease cases under the null hypothesis.
#' This is used on "Poisson" model.
#'
#' @param nb
#' A neighbors list or an adjacency matrix.
#'
#' @param name
#' A vector of names of each area.
#'
#' @param clustersize
#' The number of maximum spatial cluster size to scan, i.e., the maximum number
#' of regions included in the detected cluster
#'
#' @param maxclusters
#' The maximum number of clusters to be detected.
#'
#' @param stattype
#' Statistic type to be used (case-insensitive).
#' \describe{
#'   \item{"ORIGINAL"}{the likelihood ratio statistic by Kulldorff and
#'   Nagarwalla (1995)}
#'   \item{"RESTRICTED"}{the restricted likelihood ratio statistic by Tango
#'   (2008), with a preset parameter \code{ralpha} for restriction}
#' }
#'
#' @param scanmethod
#' Scanning method to be used (case-insensitive).
#' \describe{
#'   \item{"FLEXIBLE"}{flexible scan statistic by Tango and Takahashi (2005)}
#'   \item{"CIRCULAR"}{circular scan statistic by Kulldorff (1997)}
#' }
#'
#' @param ralpha
#' Threshold parameter of the middle p-value for the restricted likelihood ratio
#' statistic.
#'
#' @param simcount
#' The number of Monte Carlo replications to calculate a p-value for statistical
#' test.
#'
#' @param verbose
#' Print progress messages.
#'
#' @param cores
#' The number of CPU cores used for Monte Carlo replications.
#'
#'
#' @return
#' An \code{rflexscan} object which contains analysis results and specified
#' parameters.
#'
#' @references
#'   Tango T. and Takahashi K. (2005). A flexibly shaped spatial scan
#'   statistic for detecting clusters, International Journal of Health
#'   Geographics 4:11.
#'
#'   Kulldorff M. and Nagarwalla N. (1995). Spatial disease clusters:
#'   Detection and Inference. Statistics in Medicine 14:799-810.
#'
#'   Kulldorff M. (1997). A spatial scan statistic. Communications in
#'   Statistics: Theory and Methods, 26:1481-1496.
#'
#'   Tango T. (2008). A spatial scan statistic with a restricted
#'   likelihood ratio. Japanese Journal of Biometrics 29(2):75-95.
#'
#'
#' @importFrom rflexscan runFleXScan
#' @importFrom utils capture.output txtProgressBar
#' @importFrom stats BIC glm logLik poisson
#' @importFrom foreach foreach %dopar%
#' @importFrom doSNOW registerDoSNOW
#' @importFrom parallel makeCluster stopCluster
#'
#' @export
#'
spatial.scan <- function(x, y,
                         name, observed, expected, nb,
                         clustersize=15,
                         maxclusters=25,
                         stattype="ORIGINAL",
                         scanmethod="FLEXIBLE",
                         ralpha=0.2,
                         simcount=999,
                         verbose=FALSE,
                         cores = 2) {
  call <- match.call()

  stattype <- match.arg(toupper(stattype), flexscan.stattype)
  scanmethod <- match.arg(toupper(scanmethod), flexscan.scanmethod)

  # replace space
  name <- sub(" ", "_", name)

  if (!missing(x) && !missing(y)) {
    coordinates <- cbind(x, y)
    latlon <- FALSE
  } else {
    stop("Coordinates are not properly specified.")
  }

  if (missing(observed)) {
    stop("Observed numbers of diseases are not specified.")
  }

  if (!missing(expected)) {
    case <- cbind(observed, expected)
    model <- "POISSON"
  } else {
    stop("Expected numbers of diseases are not specified.")
  }

  row.names(coordinates) <- as.character(name)
  row.names(case) <- as.character(name)

  if (missing(nb)) {
    stop("A neighbours list or an adjacency matrix are not specified.")
  }

  if (is.matrix(nb)) {
    adj_mat <- nb
  } else {
    adj_mat <- matrix(0, nrow = nrow(coordinates), ncol = nrow(coordinates))
    for (i in 1:nrow(coordinates)) {
      adj_mat[i, nb[[i]]] <- 1
    }
  }
  row.names(adj_mat) <- row.names(coordinates)
  colnames(adj_mat) <- row.names(coordinates)
  diag(adj_mat) <- 2

  setting <- list()
  setting$clustersize <- clustersize
  setting$radius <- 6370
  setting$model <- as.integer(model == "BINOMIAL")
  setting$stattype <- as.integer(stattype == "RESTRICTED")
  setting$scanmethod <- as.integer(scanmethod == "CIRCULAR")
  setting$ralpha <- ralpha
  setting$cartesian <- as.integer(!latlon)
  setting$simcount <- simcount
  setting$rantype <- 0
  setting$secondary <- maxclusters - 1

  if (!verbose) {
    output <- capture.output({
      start <- date()
      clst <- runFleXScan(setting, case, coordinates, adj_mat)
      end <- date()
    })
  } else {
    start <- date()
    clst <- runFleXScan(setting, case, coordinates, adj_mat)
    end <- date()
  }

  # indicator variable z_ki を作成
  Z <- sapply(clst, function(clstr) {
    as.numeric(name %in% clstr$name)
  })
  cas_tmp <- cbind(case, as.data.frame(Z))

  neg2logLik <- numeric()
  aic <- numeric()
  bic <- numeric()
  C <- numeric()
  for (K in 0:maxclusters) {
    # Poisson regression
    retval <- glm(observed ~ . - expected, offset = log(expected), family = poisson(link = "log"), data = cas_tmp[,1:(2+K)])
    neg2logLik <- c(neg2logLik, -2 * logLik(retval))
    aic <- c(aic, retval$aic)
    bic <- c(bic, BIC(retval))
    C <- c(C, -2 * logLik(retval) + (3 * K + 1) * log(nrow(case)))
  }
  RDC <- (C[1] - C[-1]) / C[1]

  nclust <- which(RDC == max(RDC))

  # Monte Carlo test

  cl <- makeCluster(cores, type = "SOCK")
  registerDoSNOW(cl)

  pb <- txtProgressBar(max = simcount, style = 3)
  progress <- function(n) setTxtProgressBar(pb, n)
  opts <- list(progress = progress)

  setting$simcount <- 1
  maxRDC_null <- foreach (i = 1:simcount, .combine = c, .inorder = F, .options.snow = opts) %dopar% {
    sim_Observed <- sapply(case[,"expected"], function(x){
      rpois(1,x)
    })

    sim_case <- case
    sim_case[,"observed"] <- sim_Observed

    # flexscanの実行
    capture.output({
      clst_null <- runFleXScan(setting, sim_case, coordinates, adj_mat)
    })

    # indicator variable z_ki を作成
    Z <- sapply(clst_null, function(clstr) {
      as.numeric(name %in% clstr$name)
    })
    cas_tmp <- cbind(sim_case, as.data.frame(Z))

    # 情報量基準を各Kについて計算
    C_null <- numeric()
    for (K in 0:length(clst_null)) {
      # Poisson回帰
      retval <- glm(observed ~ . - expected, offset = log(expected), family = poisson(link = "log"), data = cas_tmp[,1:(2+K)])
      C_null <- c(C_null, -2 * logLik(retval) + (3 * K + 1) * log(nrow(case)))
    }
    RDC_null <- (C_null[1] - C_null[-1]) / C_null[1]

    max(RDC_null)
  }

  stopCluster(cl)

  rank <- (sum(maxRDC_null > RDC[nclust]) + 1)
  pval <- rank / (simcount + 1)

  if (toupper(model) == "POISSON") {
    colnames(case) <- c("Observed", "Expected")
  } else {
    colnames(case) <- c("Observed", "Population")
  }

  diag(adj_mat) <- 0
  for (i in 1:length(clst)) {
    x <- clst[[i]]
    adj_mat[x$area,x$area] <- adj_mat[x$area,x$area] * (10 * i)
  }

  setting$simcount <- simcount
  setting$model <- model
  setting$stattype <- stattype
  setting$scanmethod <- scanmethod
  setting$cartesian <- !latlon

  input <- list()
  input$coordinates <- coordinates
  input$case <- case
  input$adj_mat <- adj_mat

  retval <- list(call = call, input = input, cluster = clst,
                 RDC = RDC, neg2logLik = neg2logLik, AIC = aic, BIC = bic, C = C,
                 nclust = nclust, P = pval, rank = rank, null_dist = maxRDC_null,
                 setting = setting)
  class(retval) <- "multicluster"

  return(retval)
}


#' Summarizing multicluster results
#'
#' Summary method for multicluster objects.
#'
#' @param object
#' An rflexscan object to be summarized.
#'
#' @param ...
#' Ignored.
#'
#' @seealso \link{spatial.scan}
#'
#' @method summary multicluster
#' @export
#'
summary.multicluster <- function(object, ...) {
  n_cluster <- length(object$cluster)
  total_areas <- nrow(object$input$case)
  total_cases <- sum(object$input$case[,"Observed"])

  n_area <- sapply(object$cluster, function(i){length(i$area)})
  max_dist <- sapply(object$cluster, function(i) {i$max_dist})
  n_case <- sapply(object$cluster, function(i) {i$n_case})
  stats <- sapply(object$cluster, function(i) {i$stats})
  pval <- sapply(object$cluster, function(i) {i$pval})

  if (toupper(object$setting$model) == "POISSON") {
    expected <- sapply(object$cluster, function(i) {i$expected})
    RR <- sapply(object$cluster, function(i) {i$RR})

    table <- data.frame(NumArea=n_area, MaxDist=max_dist, Case=n_case,
                        Expected=expected, RR=RR, Stats=stats, P=pval)
  } else {
    population <- sapply(object$cluster, function(i) {i$population})

    table <- data.frame(NumArea=n_area, MaxDist=max_dist, Case=n_case,
                        Population=population, Stats=stats, P=pval)
  }
  row.names(table) <- 1:n_cluster

  retval <- list(call=object$call,
                 total_areas=total_areas, total_cases=total_cases,
                 cluster=table, RDC = object$RDC, nclust = object$nclust, P = object$P,
                 nclust = object$nclust, setting=object$setting)

  class(retval) <- "summary.multicluster"
  return(retval)
}


#' Print summary of multicluster results
#'
#' Print summary of multicluster results to the terminal.
#'
#' @param x
#' An summary.multicluster object to be printed.
#'
#' @param ...
#' Ignored.
#'
#' @seealso \link{spatial.scan}, \link{summary.multicluster}
#'
#' @method print summary.multicluster
#' @export
#'
print.summary.multicluster <- function(x, ...) {
  cat("\nCall:\n", paste(deparse(x$call), sep = "\n", collapse = "\n"),
      "\n\n", sep = "")

  cat("Clusters:\n")
  signif <- symnum(x$P, corr = FALSE,
                   na = FALSE, cutpoints = c(0, 0.001, 0.01, 0.05, 0.1, 1),
                   symbols = c("***", "**", "*", ".", " "))

  dig <- ceiling(log10(x$setting$simcount))

  Pm <- rep("", nrow(x$cluster))
  Pm[x$nclust] <- format(round(x$P, dig), nsmall = dig)
  sig <- rep(" ", nrow(x$cluster))
  sig[x$nclust] <- signif

  if (toupper(x$setting$model) == "POISSON") {
    table <- data.frame(NumArea = x$cluster$NumArea,
                   MaxDist = round(x$cluster$MaxDist, 3),
                   Case = x$cluster$Case,
                   Expected = round(x$cluster$Expected, 3),
                   RR = round(x$cluster$RR, 3),
                   Stats = round(x$cluster$Stats, 3),
                   Ps = format(round(x$cluster$P, dig), nsmall = dig),
                   Pm = Pm,
                   sig)
  } else {
    table <- data.frame(NumArea = x$cluster$NumArea,
                   MaxDist = round(x$cluster$MaxDist, 3),
                   Case = x$cluster$Case,
                   Population = x$cluster$Population,
                   Stats = round(x$cluster$Stats, 3),
                   Ps = format(round(x$cluster$P, dig), nsmall = dig),
                   Pm = Pm,
                   sig)
  }
  colnames(table)[ncol(table)] <- ""
  out <- capture.output(print(table, quote = FALSE, right = TRUE, print.gap = 2))
  cat(out[1:(x$nclust+1)], sep = "\n")
  cat(paste0(rep("-", nchar(out[1])), collapse = ""), "\n")
  cat(out[(x$nclust+2):length(out)], sep = "\n")
  cat("---\nSignif. codes: ", attr(signif, "legend"), "\n\n")

  cat("Number of clusters selected:", x$nclust, "\n")
  cat("")
  cat("Limit length of cluster:", x$setting$clustersize, "\n")
  cat("Number of areas:", x$total_areas, "\n")
  cat("Total cases:", x$total_cases, "\n")
  if (x$setting$cartesian) {
    cat("Coordinates: Cartesian\n")
  } else {
    cat("Coordinates: Latitude/Longitude\n")
  }
  cat("Model:", x$setting$model, "\n")
  cat("Scanning method:", x$setting$scanmethod, "\n")
  cat("Statistic type:", x$setting$stattype, "\n\n")
}
