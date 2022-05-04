#' @name threshold.selscore
#' @title Selection probabilities of the regularization model using the elastic-net penalty with multiresponse gaussian family (MNET).
#' @description
#' The penalized regression model using elastic-net penalty with multiresponse gaussian family (MNET) aims to solve the following problem:
#' \deqn{ \frac{1}{2} \| Y - X \beta \|_F^2 +  \lambda [ \frac{1-\alpha}{2} ||\beta||_F^2 + \alpha \sum_{j=1}^p ||\beta_j||_2 ]}
#'
#' @param params The parameters used in \code{main} function. If you already run the \code{main} function, it then provides the list of params.
#' @param nperm The number of permutation replicates. Default is 100.
#' @param seed An integer value for reproducible outcome. Default is 1 that does not have special meaning.
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
#' selscore <- SelectionScore(x, y, Family, seq.alpha = 0.1, seq.lambda = lambda.vec.enet, delta = 1, K = 100, psub = 1.0, penalty.factor = penalty.factor, verbose = TRUE)
#'
#' # selscore$sp$TotalScore$score %>% apply(1, max) %>% head
#' # selscore$sp$IndScore %>% apply(c(1,3), max) %>% head
#'
#' # params = selscore$params
#' # nperm = 100
#' # seed = 1
#'
#' #
#'
#'
#'
#' set.seed(1)
#' res.threshold <- threshold.selscore(selscore$params, nperm=20, perm.type = "both")
#' res.PostHoc <- selscore.PostHoc( selscore, res.threshold )
#'
#' library(dplyr)
#' FD = 10
#' with( res.PostHoc, selected[,,FD] %>% head )
#' with( res.PostHoc, selected[,,FD][selected[,"comp.total",FD]==1,] %>% head )
#' with( res.PostHoc, selected[,,FD][selected[,"sep.total",FD]==1,] %>% head )
#' with( res.PostHoc, res.PostHoc$TotalScore )
#'
#' df.PostHoc.SepTotal <- with( res.PostHoc, cbind(selected[,,FD], TotalScore=res.PostHoc$TotalScore)[selected[,"sep.total",FD]==1,] )
#' df.PostHoc.CompTotal <- with( res.PostHoc, cbind(selected[,,FD], TotalScore=res.PostHoc$TotalScore)[selected[,"comp.total",FD]==1,] )
#'
#' df.PostHoc.SepTotal[ order(df.PostHoc.SepTotal[,"TotalScore"], decreasing=TRUE), ]
#' df.PostHoc.CompTotal[ order(df.PostHoc.CompTotal[,"TotalScore"], decreasing=TRUE), ]
#'
#'
#'
#' @export threshold.selscore
threshold.selscore <- function(params, nperm = 100, seed = 1, perm.type = c("row-wise", "col-wise")) {

  # params <- selscore$params
  # params <- list(
  #         x = x,
  #         y = y,
  #         family = family,
  #         method = method,
  #         seq.alpha = seq.alpha,
  #         seq.lambda = seq.lambda,
  #         seq.df = seq.df,
  #         K = K,
  #         psub = psub,
  #         penalty.factor = penalty.factor,
  #         standardize.response = standardize.response,
  #         ...
  # )

  if( params$verbose ){
    params$verbose = FALSE
    is.verbose.converted = TRUE
  }

  n <- nrow( params$x )
  p <- sum( params$penalty.factor == 1 )
  q <- ncol( params$y )

  n.alpha <- length( params$seq.alpha )
  n.lambda <- length( params$seq.lambda )

  # empirical.result <- as.list(1:nperm)

  Empirical_Threshold <- array(0, dim=c(p, q+1) )
  maxFD = floor(0.1*p)
  Score.perm <- array(0, dim=c(maxFD, q+1, nperm, 2),
                      dimnames = list(paste0("top", 1:maxFD),
                                      c(paste0("comp.y", 1:q), "comp.total"),
                                      paste0("perm", 1:nperm),
                                      c("sort", "order")) )

  print("Permutation get started!")
  perm.params <- params
  pb <- txtProgressBar(min=0, max=nperm, style=3)
  for( perm.i in 1:nperm ){
    # set.seed((nperm+1)*j+perm.i)
    setTxtProgressBar(pb, perm.i)

    if( perm.i == 1 ) start <- proc.time()
    # permutation -------------------------------------------------------------


    # lambda.vec.enet <- grid.lambda(x = params$x, y = params$y,
    #                                family = params$family,
    #                                method = "enet",
    #                                iter = 10,
    #                                seq.alpha = params$seq.alpha, n.lambda = 10,
    #                                penalty.factor=params$penalty.factor) %>%
    #   lapply( function(x) seq(median(x), max(x), length.out=10) )


    set.seed( seed + perm.i - 1 )

    if( perm.type == "row-wise" ){
      perm.params$y <- params$y[sample(n),]

      fit.perm <- do.call("SelectionScore", perm.params)

      sc.ind <- apply( fit.perm$sp$IndScore, c(1,3), max )
      sc.total <- apply( fit.perm$sp$TotalScore, 1, max )
      sc.ind_total <- cbind(sc.ind, sc.total)
      colnames(sc.ind_total) <- c(paste0("comp.y", 1:ncol(params$y)), "comp.total")
    } else if( perm.type == "col-wise" ){
      perm.params$y <- apply(params$y, 2, function(yk) yk[sample(n)])

      fit.perm <- do.call("SelectionScore", perm.params)

      sc.ind <- apply( fit.perm$sp$IndScore, c(1,3), max )
      sc.total <- apply( fit.perm$sp$TotalScore, 1, max )
      sc.ind_total <- cbind(sc.ind, sc.total)
      colnames(sc.ind_total) <- c(paste0("sep.y", 1:ncol(params$y)), "sep.total")
    }


    # selscore$sp$IndScore %>% apply(c(1,3), max) %>% head

    # selscore$sp$TotalScore %>% apply(1, max) %>% sort %>% tail(10)
    # selscore$sp$TotalScore %>% apply(1, max) %>% order %>% tail(10)
    # fit.perm$sp$TotalScore %>% apply(1, max) %>% sort %>% tail

    # selscore$sp$IndScore %>% apply(c(1,3), max) %>% apply(2, sort, decreasing=T) %>% head
    # fit.perm$sp$IndScore %>% apply(c(1,3), max) %>% apply(2, sort, decreasing=T) %>% head

    Score.perm[,,perm.i,1] <- apply( sc.ind_total, 2, function(sck) sort(sck, decreasing=TRUE)[1:maxFD] )
    Score.perm[,,perm.i,2] <- apply( sc.ind_total, 2, function(sck) order(sck, decreasing=TRUE)[1:maxFD] )

    for( fdfd in 1:p ){
      Empirical_Threshold[fdfd, ] <-
        Empirical_Threshold[fdfd, ] +
        apply( sc.ind_total, 2, function(sck) sort(sck, decreasing=TRUE)[ fdfd ] ) / nperm
    }

    if( perm.i == 1 ) cat( "The expected time left = ", (proc.time() - start)["elapsed"] * nperm / 60, " (min) \n" )
  }

  colnames(Empirical_Threshold) <- colnames(sc.ind_total)



  return( list( threshold = Empirical_Threshold, validation = Score.perm, seed = (seed + 1:nperm - 1) ) )

}








selscore.PostHoc <- function(SelScore, Empirical_Threshold){

  if( all( names(Empirical_Threshold) == c("threshold", "validation") ) ){
    Empirical_Threshold <- Empirical_Threshold$threshold
  }

  SelScore$sp$IndScore %>% apply(c(1,3), max) %>% apply(2, sort, decreasing=TRUE) %>% head(10)
  SelScore$sp$TotalScore %>% apply(c(1), max) %>% sort(decreasing=TRUE) %>% head(10)
  Empirical_Threshold %>% head

  sel.mat <- array(NA, dim=c(p, ncol(Empirical_Threshold), p))
  dimnames(sel.mat) <- list(1:p, colnames(Empirical_Threshold), paste0("FD=", 1:p))
  for( k in 1:q ){
    for( fdfd in 1:p ){
      sel.mat[,k,fdfd] <-
        ifelse( apply(selscore$sp$IndScore, c(1,3), max)[,k] > Empirical_Threshold[fdfd,k], 1, 0 )
    }
  }

  for( k in (q+1):ncol(Empirical_Threshold) ){
    for( fdfd in 1:p ){
      sel.mat[,k,fdfd] <-
        ifelse( apply(SelScore$sp$TotalScore, 1, max) > Empirical_Threshold[fdfd,k], 1, 0 )
    }
  }





  # apply(sel.mat[,,1], 2, function(ss) sum( which(ss==1) %in% (1:p)[!(1:p) %in% 1:10] ) )
  # apply(sel.mat[,,5], 2, function(ss) sum( which(ss==1) %in% (1:p)[!(1:p) %in% 1:10] ) )
  # apply(sel.mat[,,10], 2, function(ss) sum( which(ss==1) %in% (1:p)[!(1:p) %in% 1:10] ) )
  # apply(sel.mat[,,20], 2, function(ss) sum( which(ss==1) %in% (1:p)[!(1:p) %in% 1:10] ) )
  # apply(sel.mat[,,50], 2, function(ss) sum( which(ss==1) %in% (1:p)[!(1:p) %in% 1:10] ) )

  list( selected = sel.mat, TotalScore = apply( selscore$sp$TotalScore, 1, max ) )

}


# library(mnormt)
# library(glmnet)
# library(dplyr)
#
# for( i in list.files("./R") %>% {.[!.%in%"unet-package.R"]} ){
#   source(paste0("./R/", i))
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
#                                family = Family,
#                                method = "enet",
#                                iter = 10,
#                                seq.alpha = 0.1, n.lambda = 10,
#                                penalty.factor=penalty.factor) %>%
#   lapply( function(x) seq(median(x), max(x), length.out=10) )
#
#
# selscore <- SelectionScore(x, y, Family, seq.alpha = 0.1, seq.lambda = lambda.vec.enet, delta = 1, K = 100, compare = TRUE, psub = 1.0, penalty.factor = penalty.factor, verbose = TRUE)
#
# # selscore$sp$TotalScore %>% apply(1, max) %>% head
# selscore$sp$TotalScore.withoutTune %>% apply(1, max) %>% head(10)
# selscore$sp$TotalScore.withoutWeight %>% apply(1, max) %>% head(10)
# # selscore$sp$IndScore %>% apply(c(1,3), max) %>% head
#
# selscore$sp$TotalScore %>% apply(1, max) %>% order(decreasing=TRUE) %>% head(20)
# selscore$sp$TotalScore.withoutTune %>% apply(1, max) %>% order(decreasing=TRUE) %>% head(20)
# selscore$sp$TotalScore.withoutWeight %>% apply(1, max) %>% order(decreasing=TRUE) %>% head(20)
#
#
#
#
# selscore <- SelectionScore(x, y, Family, seq.alpha = 0.1, seq.lambda = lambda.vec.enet, delta = 0.2, K = 100, compare = TRUE, psub = 1.0, penalty.factor = penalty.factor, verbose = TRUE)
#
# selscore$sp$TotalScore %>% apply(1, max) %>% order(decreasing=TRUE) %>% head(20)
# selscore$sp$TotalScore.withoutTune %>% apply(1, max) %>% order(decreasing=TRUE) %>% head(20)
# selscore$sp$TotalScore.withoutWeight %>% apply(1, max) %>% order(decreasing=TRUE) %>% head(20)
#
#
#
# # params = selscore$params
# # nperm = 100
# # seed = 1
#
# #
#
#
#
# set.seed(1)
# res.threshold <- threshold.selscore(selscore$params, nperm=5, perm.type = "row-wise")
#
#
# res.PostHoc <- selscore.PostHoc( selscore, res.threshold )
#
# library(dplyr)
# FD = 10
# with( res.PostHoc, selected[,,FD] %>% head )
# with( res.PostHoc, selected[,,FD][selected[,"comp.total",FD]==1,] %>% head )
# with( res.PostHoc, res.PostHoc$TotalScore %>% head )
#
# df.PostHoc.CompTotal <- with( res.PostHoc, cbind(selected[,,FD], TotalScore=res.PostHoc$TotalScore)[selected[,"comp.total",FD]==1,] )
#
# df.PostHoc.CompTotal[ order(df.PostHoc.CompTotal[,"TotalScore"], decreasing=TRUE), ]
