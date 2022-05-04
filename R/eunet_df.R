#' @name eunet_df
#' @title Selection probabilities of the regularization model using the elastic-net penalty with multiresponse gaussian family (MNET).
#' @description
#' The penalized regression model using elastic-net penalty with multiresponse gaussian family (MNET) aims to solve the following problem:
#' \deqn{ \frac{1}{2} \| Y - X \beta \|_F^2 +  \lambda [ \frac{1-\alpha}{2} ||\beta||_F^2 + \alpha \sum_{j=1}^p ||\beta_j||_2 ]}
#'
#' @param x Either R object or file path can be considered. A genotype data is not a data.frame but a matrix with dimension \code{p} by \code{(n+11)}. It is formatted by hapmap which has (rs, allele, chr, pos) in the first four(1-4) columns, (strand, assembly, center, protLSID, assayLSID, panel, Qcode) in the following seven(5-11) columns. If NULL, user can choose a path in interactive use.
#' @param y Either R object or file path can be considered. A phenotype data is an \code{n} by \code{p} matrix. Since the first some columns can display attributes of the phenotypes, you should enter the arguments, y.col and y.id.col, which represent the columns of phenotypes to be analyzed and the column of sample ID. If NULL, user can choose a path in interactive use.
#' @param seq.alpha Default is "object". If \code{input.type} is "object", obejects of genotype/phenotype will be entered, and if "path", paths of genotype/phenotype will be enterd. If you want to use an object, you have to make sure that the class of each column of genotype data is equal to "character".
#' @param seq.lambda The columns of phenotypes. At most 4 phenotypes can be considered, because the plot of them will be fine. Default is 2.
#' @param K The column of sample ID in the phenotype data file. Default is 1.
#' @param psub If TRUE. phenotypes are converted to be normal-shape using box-cox transformation when all phenotypes are positive.
#' @param standardize.response If TRUE. phenotypes are converted to be normal-shape using box-cox transformation when all phenotypes are positive.
#' @param penalty.factor If TRUE. phenotypes are converted to be normal-shape using box-cox transformation when all phenotypes are positive.
#'
#' @author Kipoong Kim <kkp7700@gmail.com>
#' @references
#' Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.
#' Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417-473.
#' Simon, N., Friedman, J., & Hastie, T. (2013). A blockwise descent algorithm for group-penalized multiresponse and multinomial regression. arXiv preprint arXiv:1311.6529.
#' Kim, K., Koo, J., & Sun, H. (2020). An empirical threshold of selection probability for analysis of high-dimensional correlated data. Journal of Statistical Computation and Simulation, 1-12.
#'
#' @examples
#' Data <- simdata(N=100, P=1000, Q=8, gamma = 0.5,
#'                 Y.type="continuous", n.multiclass=6,
#'                 snp.rho = 0.95, maf.min = 0.05,
#'                 y.rho = runif(1, 0.0, 0.2), scenario = "scenarioA" )
#'
#' Family <- Data$family
#'
#' lambda.list <- grid.lambda(x = Data$snp, y = Data$y,
#'                            family=Family,
#'                            iter = 10,
#'                            seq.alpha = 0.1, n.lambda = 10)
#'
#' lambda.vec <- lapply( Lambda.enet.list, function(x) seq(median(x), max(x), length.out=10) )
#' sp.mnet <- mnet(x = Data$snp, y = Data$y,
#'                 seq.alpha = 0.1,
#'                 seq.lambda = lambda.vec, K = 100)
#'
#' sp.unet <- unet(x = Data$snp, y = Data$y,
#'                 family = Family,
#'                 seq.alpha = 0.1,
#'                 seq.lambda = lambda.vec, K = 100,
#'                 penalty.factor = c(0, rep(1, ncol(Data$snp)-1)) )
#'
#' sp.unet_df <- unet_df( x = Data$snp, y = Data$y,
#'                        family = Family,
#'                        seq.alpha = 0.1,
#'                        seq.df = floor(seq(1, nrow(Data$y), length.out=10)),
#'                        K = 100, penalty.factor = c(0, rep(1, ncol(Data$snp)-1)) )
#'
#' @return A data.frame with each column representing the selection probability for each phenotype
#' @import mnormt
#' @import glmnet
#' @export eunet_df
eunet_df <-
  function(x,
           y,
           family,
           seq.alpha = NULL,
           seq.df = NULL,
           K = 100,
           psub = 0.5,
           penalty.factor = NULL,
           verbose = FALSE,
           ...) {

    # x = Data$snp
    # y = Data$y
    # family = Family
    # seq.alpha = 0.1
    # seq.df = floor(seq(1, nrow(Data$y), length.out=10))
    # K = 100
    # psub <- 0.5
    # penalty.factor <- c(0, rep(1, ncol(x)-1))

    if( length(family) != ncol(y) ) stop("The length of family should be equal to ncol(y).")
    if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
    if( is.null(seq.alpha)) seq.alpha <- 1:9*0.1
    if( is.null(penalty.factor)) penalty.factor <- rep(1, ncol(x))

    x <- as.matrix(x)
    y <- as.data.frame(y, stringsAsFactors = FALSE)

    n.pf <- sum(penalty.factor==0)
    if( n.pf > 0 ){
      wh.var <- (1:ncol(x))[-(1:n.pf)]
    } else {
      wh.var <- (1:ncol(x))
    }

    n <- nrow(x);
    p <- ncol(x)-n.pf;
    nsub <- n*psub;


    n.alpha <- length(seq.alpha)
    seq.df <- sort(seq.df, decreasing = FALSE)
    n.df <- length(seq.df)

    # IF sp -------------------------------------------------------------------
    beta.array.enet <- array( 0, dim = c(p, n.alpha, n.df, ncol(y)),
                              dimnames = list(
                                paste0("", 1:p),
                                paste0("", seq.alpha),
                                paste0("", seq.df),
                                paste0("", seq_len(ncol(y)))
                              )
    )

    names(attributes(beta.array.enet)$dimnames) <- c("Variable", "Alpha", "Df", "Phenotype")

    # beta.array.unet <-
      beta.array.unet_sum <- beta.array.unet_sum_max <- beta.array.enet


    for (kk in 1:K) {
      if( kk %% 10 == 0 ) cat("Iteration :", kk, "/", K, "\n")
      start <- proc.time()

      wsub <- sample(n, nsub)
      xsub <- x[wsub, , drop = F]
      ysub <- y[wsub, , drop = F]



      for (aa in 1:length(seq.alpha)) {

        # begin colcol ------------------------------------------------------------
        unet.matrix <- NULL
        for (colcol in 1:ncol(y)) {
          FAMILY <- family[colcol]

          if (FAMILY == "multinomial") {
            tmp.y <- model.matrix(~ . - 1,
                                  data=as.data.frame(as.factor(ysub[, colcol, drop= T]) %>% droplevels()))
            type.multinomial <- "grouped"
          } else {
            tmp.y <- ysub[, colcol, drop = T]
            type.multinomial <- NULL
          }

          fit_df_to_lambda <- df2lambda( x = xsub, y = tmp.y,
                                                   alpha = seq.alpha[aa],
                                                   # nlambda = 0.2*floor(nrow(xsub)),
                                                   # nlambda = 100,
                                                   family = FAMILY,
                                                   type.multinomial = type.multinomial,
                                                   seq.df = seq.df,
                                                   penalty.factor = penalty.factor,
                                                   verbose = verbose )

          seq.lambda <- fit_df_to_lambda$lambda


#           seq.lambda <- NULL
#           for (dfdf in seq.df) {
#             seq.lambda <- c( seq.lambda,
#                              df2lambda( x = xsub, y = tmp.y,
#                                         alpha = seq.alpha[aa],
#                                         nlambda = floor(1 * nrow(xsub)),
#                                         family = FAMILY,
#                                         type.multinomial = type.multinomial,
#                                         DFDF = dfdf,
#                                         tol = 1 ) )
#           }


          fit <- glmnet( x = xsub, y = tmp.y,
                         alpha = seq.alpha[aa],
                         lambda = seq.lambda,
                         family = FAMILY,
                         type.multinomial = type.multinomial,
                         penalty.factor = penalty.factor, ... )


          fitted.beta <-
            switch(as.character(is.list(fit$beta)),
                   "TRUE" = fit$beta[[1]],
                   "FALSE" = fit$beta)
          fitted.beta <- fitted.beta[wh.var,]

          unet.matrix <- cbind(unet.matrix, as.numeric(fitted.beta != 0))

          beta.array.enet[, aa, , colcol] <-
            beta.array.enet[, aa, , colcol] + as.numeric(fitted.beta != 0) / K
        }
        # end colcol --------------------------------------------------------------

        for (colcol in 1:ncol(y)) {
          # beta.array.unet[, aa, , colcol] <- beta.array.unet[, aa, , colcol] + apply(unet.matrix, 1, function(x) ifelse(sum(x) > 0, 1, 0)) / K
          beta.array.unet_sum[, aa, , colcol] <- beta.array.unet_sum[, aa, , colcol] + apply(unet.matrix, 1, function(x) sum(x)) / K
        }

      }

      # end alpha ---------------------------------------------------------------

      end <- proc.time()

      # cat("In this iteration, the elaspsed time=", (end - start)[3], "\n")

      if(kk == 1) cat("The expected remained time is", (end - start)[3] * (K - kk), "\n")

    }

    # end REP -----------------------------------------------------------------


    GRID <- expand.grid( lapply(dim(beta.array.unet_sum)[-1], function(x) 1:x) )
    tmp.denom1 <- array(0, c(max(GRID$Var1), max(GRID$Var2), max(GRID$Var3)))
    for (JJ in 1:nrow(GRID)) {
      # if (JJ %% 10 == 0) cat("JJ is ", JJ, "\n")
      j1 <- GRID[JJ, 1]
      j2 <- GRID[JJ, 2]
      j3 <- GRID[JJ, 3]

      tmp.denom1[j1, j2, j3] <- ifelse(max(beta.array.unet_sum[, j1, j2, j3]) == 0, 1, max(beta.array.unet_sum[, j1, j2, j3]))
      beta.array.unet_sum_max[, j1, j2, j3] <- as.array(beta.array.unet_sum[, j1, j2, j3]) / tmp.denom1[j1, j2, j3]
    }

    # beta.array.unet_sum <- beta.array.unet_sum / max(beta.array.unet_sum)

    # out <- list( enet = slam::as.simple_sparse_array( beta.array.enet ),
    #              unet = slam::as.simple_sparse_array( beta.array.unet ),
    #              unet_sum = slam::as.simple_sparse_array( beta.array.unet_sum ),
    #              unet_sum_max = slam::as.simple_sparse_array( beta.array.unet_sum_max )
    # )

    out <- list( enet = apply( beta.array.enet, c(1,4), max ),
                 # unet = apply( beta.array.unet, c(1,4), max ),
                 # unet_sum = apply( beta.array.unet_sum, c(1,4), max ),
                 unet_df = apply( beta.array.unet_sum_max, c(1,4), max )
    )


    # END sp ------------------------------------------------------------------

    return(out)

  }
