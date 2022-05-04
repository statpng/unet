#' @name SelectionScore
#' @title Selection probabilities of the regularization model using the elastic-net penalty with multiresponse gaussian family (MNET).
#' @description
#' The penalized regression model using elastic-net penalty with multiresponse gaussian family (MNET) aims to solve the following problem:
#' \deqn{ \frac{1}{2} \| Y - X \beta \|_F^2 +  \lambda [ \frac{1-\alpha}{2} ||\beta||_F^2 + \alpha \sum_{j=1}^p ||\beta_j||_2 ]}
#'
#' @param x An \code{n} by \code{p} matrix of high-dimensional genomic data such as gene expression or single nucleotide polymorphism (SNP) data.
#' @param y An \code{n} by \code{q} matrix of phenotypic outcomes, each with being either continuous or categorical.
#' @param family A vector of distribution family, such as 'gaussian', 'binomial', 'multinomial', etc.
#' @param seq.alpha A grid vector of the mixing proportion used in elastic-net penalty. Default is c(0.1, ..., 0.9).
#' @param seq.lambda A grid vector of the penalty parameter used in elastic-net penalty. It can be provided through our function 'grid.lambda'.
#' @param delta Denoising parameter ranged between 0 and 1. If \code{delta} = 1, each selection result will be tuned to have the same selection counts that is equal to the minimum of selection counts over phenotypes without any denoising process. If \code{delta} < 1, the selection counts is reduced to \code{delta} * minimum of selection counts over phenotypes.
#' @param K The number of resampling replicates when calculating the selection probability. Default is 100.
#' @param psub The subsampling proportion which reduces the computational costs efficiently. Default is 1.0 with replacement. If \code{psub} < 1.0, then resampling will be performed without replacement.
#' @param penalty.factor A vector of penalty to be imposed for each column of x. The 0 value indicates no penalty which leads the variable always to be included in the model. On the other hand, the 1 value indicates that the penalty is imposed for each column as it was.
#' @param verbose If TRUE, some information for the calculation process will be printed.
#' @param standardize.response If TRUE, each phenotype is scaled to have zero mean and unit variance.
#'
#' @author Kipoong Kim <kkp7700@gmail.com>
#' @references
#' Zou, H., & Hastie, T. (2005). Regularization and variable selection via the elastic net. Journal of the royal statistical society: series B (statistical methodology), 67(2), 301-320.
#' Meinshausen, N., & Bühlmann, P. (2010). Stability selection. Journal of the Royal Statistical Society: Series B (Statistical Methodology), 72(4), 417-473.
#' Simon, N., Friedman, J., & Hastie, T. (2013). A blockwise descent algorithm for group-penalized multiresponse and multinomial regression. arXiv preprint arXiv:1311.6529.
#' Kim, K., Koo, J., & Sun, H. (2020). An empirical threshold of selection probability for analysis of high-dimensional correlated data. Journal of Statistical Computation and Simulation, 1-12.
#'
#'
#' @return A data.frame with each column representing the selection probability for each phenotype
#' @import mnormt
#' @import glmnet
#'
#' @examples
#' # library(mnormt)
#' library(mnormt)
#' library(glmnet)
#' library(dplyr)
#'
#' for( i in list.files("./R") %>% {.[!.%in%"unet-package.R"]} ){
#'   source(paste0("./R/", i))
#' }
#'
#' set.seed(1)
#' n=100; p=1000; q=4
#'
#' X <- replicate(p, rbinom(n,2,0.2)) #generate the SNP X=0,1,2
#' b <- matrix(0, p, q)
#' b[1:5,1:4] <- 1.0
#' b[6:10,1:2] <- 1.0
#'
#'
#' Z <- replicate(1, rnorm(n))
#' g <- matrix(0, 1, q)
#' g[1,1:4] <- 0.1
#'
#' x <- cbind(Z, X)
#' beta <- rbind(g, b)
#'
#' y <- x%*%beta + replicate(q, rnorm(n,1))
#' y[,2:3] <- apply(y[,2:3], 2, function(yk) ifelse(yk > median(yk), 1, 0) )
#' Family<-c("gaussian","binomial","binomial","gaussian")
#' alpha.vec <- 0.1
#' penalty.factor <- c(0, rep(1, p))
#'
#' lambda.vec.enet <- grid.lambda(x = x, y = y,
#'                                family = Family,
#'                                method = "enet",
#'                                iter = 10,
#'                                seq.alpha = 0.1, n.lambda = 10,
#'                                penalty.factor=penalty.factor) %>%
#'   lapply( function(x) seq(median(x), max(x), length.out=10) )
#'
#'
#' SelectionScore(x, y, Family, seq.alpha = 0.1, seq.lambda = lambda.vec.enet, delta = 1, K = 100, psub = 1.0, penalty.factor = penalty.factor, verbose = TRUE)
#'
#'
#'
#' set.seed(1)
#' sp.total2 <- main( x = cbind(Z,X), y = Y,
#' family = Family,
#' method = c("unet_df"),
#' seq.alpha = alpha.vec,
#' seq.lambda = list(enet=NULL,
#' mnet=NULL),
#' seq.df = df.vec,
#' K = 10,
#' penalty.factor = penalty.factor )
#'
#' set.seed(1)
#' res.threshold2 <- threshold(sp.total2$params, nperm=10)
#'
#'
#'
#'
#'
#' @export SelectionScore
SelectionScore <- function(x,
                         y,
                         family,
                         seq.alpha = NULL,
                         seq.lambda = NULL,
                         delta = 1.0,
                         compare = TRUE,
                         K = 100,
                         psub = 1.0,
                         penalty.factor = NULL,
                         standardize.response = TRUE,
                         verbose = FALSE,
                         # const = NULL,
                         ...) {

  # seq.alpha = NULL
  # seq.lambda = NULL
  # delta = 1.0
  # compare = TRUE
  # K = 100
  # psub = 1.0
  # penalty.factor = NULL
  # standardize.response = TRUE
  # verbose = FALSE
  # ... = NULL


  # library(mnormt)
  # library(glmnet)
  # library(dplyr)
  # for( i in list.files("./R") %>% {.[!.%in%"unet-package.R"]} ){
  # source(paste0("./R/", i))
  # }
  #
  # set.seed(1)
  # n=100; p=1000; q=4
  #
  # X <- replicate(p, rbinom(n,2,0.2)) #generate the SNP X=0,1,2
  # b <- matrix(0, p, q)
  # b[1:5,1:4] <- 1.0
  # b[6:10,1:2] <- 1.0
  #
  #
  # Z <- replicate(1, rnorm(n))
  # g <- matrix(0, 1, q)
  # g[1,1:4] <- 0.1
  #
  # x <- cbind(Z, X)
  # beta <- rbind(g, b)
  #
  # y <- x%*%beta + replicate(q, rnorm(n,1))
  # y[,2:3] <- apply(y[,2:3], 2, function(yk) ifelse(yk > median(yk), 1, 0) )
  # Family<-c("gaussian","binomial","binomial","gaussian")
  # alpha.vec <- 0.1
  # penalty.factor <- c(0, rep(1, p))
  #
  # lambda.vec.enet <- grid.lambda(x = x, y = y,
  #                                 family = Family,
  #                                 method = "enet",
  #                                 iter = 10,
  #                                 seq.alpha = 0.1, n.lambda = 10,
  #                                 penalty.factor=penalty.factor) %>%
  #   lapply( function(x) seq(median(x), max(x), length.out=10) )
  #
  #
  # x = x
  # y = y
  # family = Family
  # seq.alpha = alpha.vec
  # seq.lambda = lambda.vec.enet
  # K = 100
  # psub = 1.0
  # penalty.factor = penalty.factor
  # standardize.response = TRUE
  # verbose = FALSE
  # const = NULL
  # ...=NULL



  # params <- perm.params
  # x <- params$x
  # y <- params$y
  # family <- params$family
  # seq.alpha <- params$seq.alpha
  # seq.lambda <- params$seq.lambda
  # delta <- params$delta
  # K <- params$K
  # psub <- params$psub
  # penalty.factor <- params$penalty.factor
  # standardize.response <- params$standardize.response
  # verbose <- params$verbose



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

  sample.replace = FALSE
  if( psub == 1.0 ) sample.replace = TRUE

  n.alpha <- length(seq.alpha)
  seq.lambda <- lapply( seq.lambda, sort, decreasing = FALSE )
  n.lambda <- min( sapply(seq.lambda, length) )



  beta.array.enet <- array( 0, dim = c(p, n.alpha, n.lambda, ncol(y)),
        dimnames = list(
                paste0("", 1:p),
                paste0("", seq.alpha),
                paste0("", seq_len(n.lambda)),
                paste0("", seq_len(ncol(y)))
        ) )
  names(attributes(beta.array.enet)$dimnames) <- c("Variable", "Alpha", "Lambda", "Phenotype")




  gc()


  for (kk in 1:K) {
    if( kk %% 10 == 0 & verbose ) cat("Iteration :", kk, "/", K, "\n")
    start <- proc.time()

    wsub <- sample(n, nsub, replace = sample.replace)
    xsub <- x[wsub, , drop = F]
    ysub <- y[wsub, , drop = F]


    for (aa in 1:length(seq.alpha)) {

      for (colcol in 1:ncol(y)) {

        FAMILY <- family[colcol]

        type.multinomial <- NULL
        tmp.y <- ysub[, colcol, drop= T]

        if (FAMILY == "multinomial"){
          tmp.y <- y.multinom(ysub[, colcol, drop= T])
          type.multinomial <- "grouped"
        }


        fit <- glmnet( x = xsub, y = tmp.y,
                       alpha = seq.alpha[aa],
                       lambda = seq.lambda[[colcol]],
                       family = FAMILY,
                       type.multinomial = type.multinomial,
                       standardize.response = standardize.response,
                       penalty.factor = penalty.factor, ... )

        fitted.beta.enet <-
          switch(as.character(is.list(fit$beta)),
                 "TRUE" = fit$beta[[1]],
                 "FALSE" = fit$beta)
        fitted.beta.enet <- fitted.beta.enet[wh.var,]


        beta.array.enet[, aa, , colcol] <- beta.array.enet[, aa, , colcol] + as.numeric(fitted.beta.enet != 0)

      }

    }

    end <- proc.time()

    if(kk == 1 & verbose) cat("\n The expected remained time is", (end - start)[3] * (K - kk), "\n")

  }

  # beta.array.enet[,1,10,] %>% head(20)
  # beta.array.enet[,1,10,1] %>% table
  # hist(beta.array.enet[,1,10,1])



  SCORE <- array(NA, dim(beta.array.enet)[-2])
  SCORE.withoutWeight <- array(NA, dim(beta.array.enet)[-2])
  SCORE.withoutTune <- array(NA, dim(beta.array.enet)[-2])

  SelScore <- SelScore.withoutTune <- SelScore.withoutWeight <- NULL
  for( lamb in 1:n.lambda ){
    SCORE[,lamb,] <- png.enet2score( beta.array.enet[,1,lamb,], delta = delta )$score
    if( compare ){
      SCORE.withoutTune[,lamb,] <- beta.array.enet[,1,lamb,]
      SCORE.withoutWeight[,lamb,] <- png.enet2score( beta.array.enet[,1,lamb,], delta = delta, weight = FALSE )$score
    }

    SelScore <- cbind(SelScore, apply( SCORE[,lamb,], 1, sum ) / K)
    if( compare ){
      SelScore.withoutTune <- cbind(SelScore.withoutTune, apply( SCORE.withoutTune[,lamb,], 1, sum ) / K)
      SelScore.withoutWeight <- cbind(SelScore.withoutWeight, apply( SCORE.withoutWeight[,lamb,], 1, sum ) / K)
    }

  }
  colnames(SelScore) <- paste0("lambda.", 1:n.lambda)

  out <- SelScore






  # if( is.null(const) ){
  #   MaxConst <- apply( beta.array.unet_sum, 2:4, function(x){
  #     ifelse(max(x) == 0, 1, max(x))
  #   } )
  #
  #   beta.array.unet_df <- apply( beta.array.unet_sum, 2:4, function(x){
  #     MAX <- ifelse(max(x) == 0, 1, max(x))
  #     x / MAX
  #   } )
  #
  #   const <- MaxConst
  # } else {
  #
  #   GRID <- expand.grid( lapply(dim(beta.array.unet_sum)[-1], function(x) 1:x) )
  #   tmp.denom1 <- array(0, c(max(GRID$Var1), max(GRID$Var2), max(GRID$Var3)))
  #   for (JJ in 1:nrow(GRID)) {
  #     # if (JJ %% 10 == 0) cat("JJ is ", JJ, "\n")
  #     j1 <- GRID[JJ, 1]
  #     j2 <- GRID[JJ, 2]
  #     j3 <- GRID[JJ, 3]
  #
  #     beta.array.unet_df[, j1, j2, j3] <- as.array(beta.array.unet_sum[, j1, j2, j3]) / const[j1, j2, j3]
  #   }
  #
  # }




  params <- list(
    x = x,
    y = y,
    family = family,
    seq.alpha = seq.alpha,
    seq.lambda = seq.lambda,
    delta = delta,
    K = K,
    psub = psub,
    penalty.factor = penalty.factor,
    standardize.response = standardize.response,
    verbose = verbose,
    # const = const,
    ...
  )


  if( compare ){
    colnames(SelScore.withoutTune) <- colnames(SelScore.withoutWeight) <- paste0("lambda.", 1:n.lambda)

    return(list(sp = list(IndScore=SCORE/K,
                          TotalScore=out,
                          TotalScore.withoutTune = SelScore.withoutTune,
                          TotalScore.withoutWeight = SelScore.withoutWeight), params = params))

  } else {

    return(list(sp = list(IndScore=SCORE/K, TotalScore=out), params = params))

  }



}
