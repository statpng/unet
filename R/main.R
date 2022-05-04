#' @name main
#' @title Selection probabilities of the regularization model using the elastic-net penalty with multiresponse gaussian family (MNET).
#' @description
#' The penalized regression model using elastic-net penalty with multiresponse gaussian family (MNET) aims to solve the following problem:
#' \deqn{ \frac{1}{2} \| Y - X \beta \|_F^2 +  \lambda [ \frac{1-\alpha}{2} ||\beta||_F^2 + \alpha \sum_{j=1}^p ||\beta_j||_2 ]}
#'
#' @param x An \code{n} by \code{p} matrix of high-dimensional genomic data such as gene expression or single nucleotide polymorphism (SNP) data.
#' @param y An \code{n} by \code{q} matrix of phenotypic outcomes, each with being either continuous or categorical.
#' @param family A vector of distribution family, such as 'gaussian', 'binomial', 'multinomial', etc.
#' @param method A vector of methods with 'enet', 'mnet', 'unet' and 'unet_df'
#' @param seq.alpha A grid vector of the mixing proportion used in elastic-net penalty. Default is c(0.1, ..., 0.9).
#' @param seq.lambda A grid vector of the penalty parameter used in elastic-net penalty. It can be provided through our function 'grid.lambda'.
#' @param K The number of resampling replicates when calculating the selection probability. Default is 100.
#' @param psub The subsampling proportion which reduces the computational costs efficiently. Default is 0.5.
#' @param penalty.factor A vector of penalty to be imposed for each column of x. The 0 value indicates no penalty which leads the variable always to be included in the model. On the other hand, the 1 value indicates that the penalty is imposed for each column as it was.
#' @param verbose If TRUE, some information for the calculation process will be printed.
#' @param standardize.response If TRUE, each phenotype is scaled to have zero mean and unit variance.
#' @param const Normalizing Constant which makes the selection probability to be in [0, 1]. Don't touch this argument.
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
#' # library(glmnet)
#' # library(dplyr)
#' # for( i in list.files("./R") %>% {.[!.%in%"unet-package.R"]} ){
#' #        source(paste0("./R/", i))
#' # }
#'
#'
#'
#' set.seed(1)
#'
#' n=100; p=1000; q=4
#'
#' X <- replicate(p, rbinom(n,2,0.2)) #generate the SNP X=0,1,2
#' b <- matrix(0, p, q)
#' b[1:5,1:2] <- 1.0
#'
#'
#' Z <- replicate(1, rnorm(n))
#' g <- matrix(0, 1, q)
#' g[1,1:4] <- 0.1
#'
#'
#' Y <- cbind(Z, X)%*%rbind(g, b) + replicate(q, rnorm(n,1))
#' Y[,2:3] <- apply(Y[,2:3], 2, function(yk) ifelse(yk > median(yk), 1, 0) )
#' Family<-c("gaussian","binomial","binomial","gaussian")
#' alpha.vec <- 1:9*0.1
#' df.vec = floor(seq(1, floor(n/2), length.out=10))
#' penalty.factor <- c(0, rep(1, p))
#'
#'
#'
#'
#' lambda.list.enet <- grid.lambda(x = X, y = Y,
#' family = Family,
#' method = "enet",
#' iter = 10,
#' seq.alpha = 1:9*0.1, n.lambda = 10,
#' penalty.factor=penalty.factor)
#'
#'
# lambda.vec.enet <- lapply( lambda.list.enet, function(x) seq(median(x), max(x), length.out=10) )
#'
#'
#'
#'
#'
#'
#'
#' set.seed(1)
#' sp.unet <- unet(x = cbind(Z,X), y = Y,
#' family = Family,
#' seq.alpha = alpha.vec,
#' seq.lambda = lambda.vec.enet, K = 10,
#' penalty.factor = penalty.factor )
#'
#'
#' set.seed(1)
#' sp.unet_df <- unet_df( x = cbind(Z,X), y = Y,
#' family = Family,
#' seq.alpha = alpha.vec,
#' seq.df = df.vec,
#' K = 10,
#' penalty.factor = penalty.factor )
#'
#'
#' set.seed(1)
#' sp.total <- main( x = cbind(Z,X), y = Y,
#' family = Family,
#' method = c("enet", "unet", "unet_df"),
#' seq.alpha = alpha.vec,
#' seq.lambda = list(enet=lambda.vec.enet,
#' mnet=NULL),
#' seq.df = df.vec,
#' K = 10,
#' penalty.factor = penalty.factor )
#'
#' SelProb <- sp.total$sp
#'
#' head(SelProb)
#'
#'
#'
#'
#' set.seed(1)
#' sp.total1 <- main( x = cbind(Z,X), y = Y,
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
#' res.threshold1 = threshold(sp.total$params, nperm=10)
#'
#'
#'
#'
#'
#' set.seed(1)
#'
#' n=100; p=1000; q=4
#'
#' X <- replicate(p, rbinom(n,2,0.2)) #generate the SNP X=0,1,2
#' b <- matrix(0, p, q)
#'
#'
#' Z <- replicate(1, rnorm(n))
#' g <- matrix(0, 1, q)
#'
#'
#' Y <- cbind(Z, X)%*%rbind(g, b) + replicate(q, rnorm(n,1))
#' Y[,2:3] <- apply(Y[,2:3], 2, function(yk) ifelse(yk > median(yk), 1, 0) )
#' Family<-c("gaussian","binomial","binomial","gaussian")
#' alpha.vec <- 1:9*0.1
#' df.vec = floor(seq(1, floor(n/2), length.out=10))
#' penalty.factor <- c(0, rep(1, p))
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
#' @export main
main <- function(x,
                 y,
                 family,
                 method,
                 seq.alpha = NULL,
                 seq.lambda = NULL,
                 seq.df = NULL,
                 K = 100,
                 psub = 0.5,
                 penalty.factor = NULL,
                 standardize.response = TRUE,
                 verbose = FALSE,
                 const = NULL,
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
        if( is.null(seq.df)) seq.df <- floor(seq(1, floor(nrow(y)), length.out=10))
        if( is.null(penalty.factor)) penalty.factor <- rep(1, ncol(x))

        include.enet <- include.mnet <- include.unet <- include.unet_df <- FALSE
        if( "enet" %in% method ) include.enet <- TRUE
        if( "mnet" %in% method ) include.mnet <- TRUE
        if( "unet" %in% method ) include.unet <- TRUE
        if( "unet_df" %in% method ) include.unet_df <- TRUE

        include.lambda.enet <- include.lambda.mnet <- TRUE
        if( is.null(seq.lambda$enet) ) include.lambda.enet <- FALSE
        if( is.null(seq.lambda$mnet) ) include.lambda.mnet <- FALSE
        # if( (!include.lambda.enet) & (!include.lambda.mnet) ) stop("Either 'enet' or 'mnet' element should be included in the lambda.seq.")

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

        if(include.enet | include.unet | include.unet_df){
                seq.lambda.enet <- seq.lambda$enet
                seq.lambda.enet <- lapply( seq.lambda.enet, sort, decreasing = FALSE )
                n.lambda.enet <- unique( sapply(seq.lambda.enet, length) )
        }
        if(include.mnet){
                seq.lambda.mnet <- seq.lambda$mnet
                seq.lambda.mnet <- lapply( seq.lambda.mnet, sort, decreasing = FALSE )
                n.lambda.mnet <- unique( sapply(seq.lambda.mnet, length) )
        }

        if(include.unet_df){
                seq.df <- sort(seq.df, decreasing = FALSE)
                n.df <- length(seq.df)
        }




        if( include.enet ){
                beta.array.enet <- array( 0, dim = c(p, n.alpha, n.lambda.enet, ncol(y)),
                                         dimnames = list(
                                                 paste0("", 1:p),
                                                 paste0("", seq.alpha),
                                                 paste0("", seq_len(n.lambda.enet)),
                                                 paste0("", seq_len(ncol(y)))
                                         ) )
                names(attributes(beta.array.enet)$dimnames) <- c("Variable", "Alpha", "Lambda", "Phenotype")
        }
        if( include.mnet ){
                beta.array.mnet <- array( 0, dim = c(p, n.alpha, n.lambda.mnet, ncol(y)),
                                          dimnames = list(
                                                  paste0("", 1:p),
                                                  paste0("", seq.alpha),
                                                  paste0("", seq_len(n.lambda.mnet)),
                                                  paste0("", seq_len(ncol(y)))
                                          ) )
                names(attributes(beta.array.mnet)$dimnames) <- c("Variable", "Alpha", "Lambda", "Phenotype")
        }
        if( include.unet ){
                beta.array.unet <- array( 0, dim = c(p, n.alpha, n.lambda.enet, ncol(y)),
                                          dimnames = list(
                                                  paste0("", 1:p),
                                                  paste0("", seq.alpha),
                                                  paste0("", seq_len(n.lambda.enet)),
                                                  paste0("", seq_len(ncol(y)))
                                          ) )
                names(attributes(beta.array.unet)$dimnames) <- c("Variable", "Alpha", "Lambda", "Phenotype")
        }
        if( include.unet_df ){
                beta.array.unet_df <- array( 0, dim = c(p, n.alpha, n.df, ncol(y)),
                                         dimnames = list(
                                                 paste0("", 1:p),
                                                 paste0("", seq.alpha),
                                                 paste0("", seq.df),
                                                 paste0("", seq_len(ncol(y)))
                                         ) )
                names(attributes(beta.array.unet_df)$dimnames) <- c("Variable", "Alpha", "Df", "Phenotype")

                beta.array.unet_sum <- beta.array.unet_df
        }

        gc()

        for (kk in 1:K) {
                if( kk %% 10 == 0 & verbose ) cat("Iteration :", kk, "/", K, "\n")
                start <- proc.time()

                wsub <- sample(n, nsub)
                xsub <- x[wsub, , drop = F]
                ysub <- y[wsub, , drop = F]


                for (aa in 1:length(seq.alpha)) {

                        # begin colcol ------------------------------------------------------------
                        unet.matrix <- unet_df.matrix <- NULL
                        for (colcol in 1:ncol(y)) {
                                FAMILY <- family[colcol]

                                type.multinomial <- NULL
                                tmp.y <- ysub[, colcol, drop= T]

                                if (FAMILY == "multinomial"){
                                  tmp.y <- y.multinom(ysub[, colcol, drop= T])
                                  type.multinomial <- "grouped"
                                }


                                # enet --------------------------------------------------------------------

                                if(include.enet | include.unet){
                                        fit <- glmnet( x = xsub, y = tmp.y,
                                                       alpha = seq.alpha[aa],
                                                       lambda = seq.lambda.enet[[colcol]],
                                                       family = FAMILY,
                                                       type.multinomial = type.multinomial,
                                                       standardize.response = standardize.response,
                                                       penalty.factor = penalty.factor, ... )

                                        fitted.beta.enet <-
                                                switch(as.character(is.list(fit$beta)),
                                                       "TRUE" = fit$beta[[1]],
                                                       "FALSE" = fit$beta)
                                        fitted.beta.enet <- fitted.beta.enet[wh.var,]
                                }



                                if(include.enet){
                                  beta.array.enet[, aa, , colcol] <- beta.array.enet[, aa, , colcol] + as.numeric(fitted.beta.enet != 0) / K
                                }

                                if(include.unet){
                                  unet.matrix <- cbind(unet.matrix, as.numeric(fitted.beta.enet != 0))
                                }


                                if(include.unet_df){
                                  fit_df_to_lambda <- df2lambda( x = xsub, y = tmp.y,
                                                                 alpha = seq.alpha[aa],
                                                                 family = FAMILY,
                                                                 type.multinomial = type.multinomial,
                                                                 seq.df = seq.df,
                                                                 penalty.factor = penalty.factor,
                                                                 verbose = verbose )

                                  seq.dflambda <- fit_df_to_lambda$lambda

                                  fit <- glmnet( x = xsub, y = tmp.y,
                                                 alpha = seq.alpha[aa],
                                                 lambda = seq.dflambda,
                                                 family = FAMILY,
                                                 type.multinomial = type.multinomial,
                                                 standardize.response = standardize.response,
                                                 penalty.factor = penalty.factor, ... )

                                  fitted.beta <-
                                          switch(as.character(is.list(fit$beta)),
                                                 "TRUE" = fit$beta[[1]],
                                                 "FALSE" = fit$beta)
                                  fitted.beta <- fitted.beta[wh.var,]

                                  unet_df.matrix <- cbind(unet_df.matrix, as.numeric(fitted.beta != 0))
                                }


                                # mnet --------------------------------------------------------------------

                                if(include.mnet){
                                  fit <- glmnet( x = xsub, y = tmp.y,
                                                 alpha = seq.alpha[aa],
                                                 lambda = seq.lambda.mnet[[colcol]],
                                                 family = "mgaussian",
                                                 type.multinomial = type.multinomial,
                                                 standardize.response = standardize.response,
                                                 penalty.factor = penalty.factor, ... )

                                  fitted.beta <-
                                          switch(as.character(is.list(fit$beta)),
                                                 "TRUE" = fit$beta[[1]],
                                                 "FALSE" = fit$beta)
                                  fitted.beta <- fitted.beta[wh.var,]

                                  beta.array.mnet[,aa,,colcol] <- beta.array.mnet[,aa,,colcol] + as.numeric( fitted.beta != 0 )/K
                                }

                        }

                        for( colcol in 1:ncol(y) ){
                          if( include.unet ){
                            beta.array.unet[,aa,,colcol] <- beta.array.unet[,aa,,colcol] + apply( unet.matrix, 1, function(x) ifelse( sum(x)>0, 1, 0 ) )/K
                          }

                          if( include.unet_df ){
                            beta.array.unet_sum[,aa,,colcol] <- beta.array.unet_sum[,aa,,colcol] + apply(unet_df.matrix, 1, function(x) sum(x)) / K
                          }
                        }

                }

                # end alpha ---------------------------------------------------------------

                end <- proc.time()

                if(kk == 1 & verbose) cat("\n The expected remained time is", (end - start)[3] * (K - kk), "\n")

        }

        # end REP -----------------------------------------------------------------


        if( include.unet_df ){

          # GRID <- expand.grid( lapply(dim(beta.array.unet_sum)[-1], function(x) 1:x) )
          # tmp.denom1 <- array(0, c(max(GRID$Var1), max(GRID$Var2), max(GRID$Var3)))
          # for (JJ in 1:nrow(GRID)) {
          #         # if (JJ %% 10 == 0) cat("JJ is ", JJ, "\n")
          #         j1 <- GRID[JJ, 1]
          #         j2 <- GRID[JJ, 2]
          #         j3 <- GRID[JJ, 3]
          #
          #         tmp.denom1[j1, j2, j3] <- ifelse(max(beta.array.unet_sum[, j1, j2, j3]) == 0, 1, max(beta.array.unet_sum[, j1, j2, j3]))
          #         beta.array.unet_df[, j1, j2, j3] <- as.array(beta.array.unet_sum[, j1, j2, j3]) / tmp.denom1[j1, j2, j3]
          # }

                if( is.null(const) ){
                        MaxConst <- apply( beta.array.unet_sum, 2:4, function(x){
                                ifelse(max(x) == 0, 1, max(x))
                        } )

                        beta.array.unet_df <- apply( beta.array.unet_sum, 2:4, function(x){
                                MAX <- ifelse(max(x) == 0, 1, max(x))
                                x / MAX
                        } )

                        const <- MaxConst
                } else {

                        GRID <- expand.grid( lapply(dim(beta.array.unet_sum)[-1], function(x) 1:x) )
                        tmp.denom1 <- array(0, c(max(GRID$Var1), max(GRID$Var2), max(GRID$Var3)))
                        for (JJ in 1:nrow(GRID)) {
                                # if (JJ %% 10 == 0) cat("JJ is ", JJ, "\n")
                                j1 <- GRID[JJ, 1]
                                j2 <- GRID[JJ, 2]
                                j3 <- GRID[JJ, 3]

                                beta.array.unet_df[, j1, j2, j3] <- as.array(beta.array.unet_sum[, j1, j2, j3]) / const[j1, j2, j3]
                        }

                }


        }

        out <- list(NULL)
        if( include.enet ) out <- append(out, list(enet=apply( beta.array.enet, c(1,4), max) ))
        if( include.mnet ) out <- append(out, list(mnet=apply( beta.array.mnet, c(1,4), max) ))
        if( include.unet ) out <- append(out, list(unet=apply( beta.array.unet, c(1,4), max) ))
        if( include.unet_df ){
        	out <- append(out, list(unet_sumdf=apply( beta.array.unet_sum, c(1,4), max) ))
        	out <- append(out, list(unet_df=apply( beta.array.unet_df, c(1,4), max) ))
        }


        out <- Filter( Negate(is.null), out )

        # END sp ------------------------------------------------------------------

        params <- list(
                x = x,
                y = y,
                family = family,
                method = method,
                seq.alpha = seq.alpha,
                seq.lambda = seq.lambda,
                seq.df = seq.df,
                K = K,
                psub = psub,
                penalty.factor = penalty.factor,
                standardize.response = standardize.response,
                verbose = verbose,
                const = const,
                ...
        )

        return(list(sp = out, params = params))

}




# library(mnormt)
# library(glmnet)
# library(dplyr)
# for( i in list.files("./R") %>% {.[!.%in%"unet-package.R"]} ){
#        source(paste0("./R/", i))
# }
#
# set.seed(1)
#
# n=100; p=1000; q=4
#
# X <- replicate(p, rbinom(n,2,0.2)) #generate the SNP X=0,1,2
# b <- matrix(0, p, q)
# b[1:5,1:2] <- 1.0
#
#
# Z <- replicate(1, rnorm(n))
# g <- matrix(0, 1, q)
# g[1,1:4] <- 0.1
#
#
# Y <- cbind(Z, X)%*%rbind(g, b) + replicate(q, rnorm(n,1))
# Y[,2:3] <- apply(Y[,2:3], 2, function(yk) ifelse(yk > median(yk), 1, 0) )
# Family<-c("gaussian","binomial","binomial","gaussian")
# alpha.vec <- 1:9*0.1
# df.vec = floor(seq(1, floor(n/2), length.out=10))
# penalty.factor <- c(0, rep(1, p))
#
#
#
#
# lambda.list.enet <- grid.lambda(x = X, y = Y,
# family = Family,
# method = "enet",
# iter = 10,
# seq.alpha = 1:9*0.1, n.lambda = 10,
# penalty.factor=penalty.factor)
#
# lambda.vec.enet <- lapply( lambda.list.enet, function(x) seq(median(x), max(x), length.out=10) )
#
# set.seed(1)
# sp.total <- main( x = cbind(Z,X), y = Y,
#                         family = Family,
#                         method = c("unet_df"),
#                         seq.alpha = alpha.vec,
#                         seq.lambda = list(enet=lambda.vec.enet,
#                         mnet=NULL),
#                         seq.df = df.vec,
#                         K = 20,
#                         penalty.factor = penalty.factor )
#
# SelProb <- sp.total$sp$unet_df[,1]
#
# head(SelProb)
#
# set.seed(1)
# res.threshold1 = threshold(sp.total$params, type="row-wise", nperm=10)
# set.seed(1)
# res.threshold2 = threshold(sp.total$params, type="col-wise", nperm=10)
#
# res.threshold1$unet_df[,1][10]
# res.threshold2$unet_df[,1][10]
#
# sort(SelProb) %>% tail(20)
# #
#
#
# set.seed(1)
# sp.total1 <- main( x = cbind(Z,X), y = Y,
# family = Family,
# method = c("unet_df"),
# seq.alpha = alpha.vec,
# seq.lambda = list(enet=NULL,
# mnet=NULL),
# seq.df = df.vec,
# K = 10,
# penalty.factor = penalty.factor )
#
# set.seed(1)
# res.threshold1 = threshold(sp.total$params, nperm=10)
