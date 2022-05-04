#' @name unet_df
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
#' penalty.factor <- c(0, rep(1, p))
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
#' lambda.vec.enet <- lapply( lambda.list.enet, function(x) seq(median(x), max(x), length.out=10) )
#'
#'
#'
#' alpha.vec <- 1:9*0.1
#' df.vec = floor(seq(1, floor(n/2), length.out=10))
#' penalty.factor <- c(0, rep(1, p))
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
#' @export unet_df
unet_df <- function(x,
                 y,
                 family,
                 seq.alpha = NULL,
                 seq.df = NULL,
                 K = 100,
                 psub = 0.5,
                 penalty.factor = NULL,
                 standardize.response = TRUE,
                 verbose = TRUE,
                 ...) {

        # x = Data$snp
        # y = Data$y
        # family = Family
        # seq.alpha = 0.1
        # seq.df = floor(seq(1, nrow(Data$y), length.out=10))
        # K = 100
        # psub <- 0.5
        # penalty.factor <- c(0, rep(1, ncol(x)-1))

        # library(glmnet)
        # x = cbind(Z,X)
        # y = Y
        # family = Family
        # method = c("enet", "unet", "unet_df")
        # seq.alpha = 1:9*0.1
        # seq.df = NULL
        # K = 10
        # psub = 0.5
        # penalty.factor = penalty.factor
        # standardize.response = TRUE
        # verbose=TRUE


        # \begin{Data} ------------------------------------------------------------

        if( length(family) != ncol(y) ) stop("The length of family should be equal to ncol(y).")
        if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")
        if( is.null(seq.alpha)) seq.alpha <- 1:9*0.1
        if( is.null(seq.df)) seq.df <- floor(seq(1, floor(nrow(y)), length.out=10))
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

        # \end{Data} --------------------------------------------------------------



        seq.df <- sort(seq.df, decreasing = FALSE)
        n.df <- length(seq.df)



        beta.array.unet_df <- array( 0, dim = c(p, n.alpha, n.df),#, ncol(y)),
                                     dimnames = list(
                                             paste0("", 1:p),
                                             paste0("", seq.alpha),
                                             paste0("", seq.df)#,
                                             # paste0("", seq_len(ncol(y)))
                                     ) )
        names(attributes(beta.array.unet_df)$dimnames) <- c("Variable", "Alpha", "Df")#, "Phenotype")

        beta.array.unet_sum <- beta.array.unet_df



        gc()



        for (kk in 1:K) {
                if( kk %% 10 == 0 ) cat("Iteration :", kk, "/", K, "\n")
                start <- proc.time()

                wsub <- sample(n, nsub)
                xsub <- x[wsub, , drop = F]
                ysub <- y[wsub, , drop = F]


                for (aa in 1:length(seq.alpha)) {


                        # begin colcol ------------------------------------------------------------
                        unet_df.matrix <- NULL
                        for (colcol in 1:ncol(y)) {
                                FAMILY <- family[colcol]

                                type.multinomial <- NULL
                                tmp.y <- ysub[, colcol, drop= T]

                                if (FAMILY == "multinomial"){
                                        tmp.y <- y.multinom(ysub[, colcol, drop= T])
                                        type.multinomial <- "grouped"
                                }


                                # enet --------------------------------------------------------------------


                                fit_df_to_lambda <- df2lambda( x = xsub, y = tmp.y,
                                                               alpha = seq.alpha[aa],
                                                               family = FAMILY,
                                                               type.multinomial = type.multinomial,
                                                               seq.df = seq.df,
                                                               penalty.factor = penalty.factor,
                                                               verbose = FALSE )

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

                        beta.array.unet_sum[,aa,] <- beta.array.unet_sum[,aa,] + apply(unet_df.matrix, 1, function(x) sum(x)) / K

                }

                # end alpha ---------------------------------------------------------------

                end <- proc.time()

                if(kk == 1 & verbose) cat("\n The expected remained time is", (end - start)[3] * (K - kk), "\n")

        }

        # end REP -----------------------------------------------------------------


        MaxConst <- apply( beta.array.unet_sum, 2:3, function(x){
                ifelse(max(x) == 0, 1, max(x))
        } )

        beta.array.unet_df <- apply( beta.array.unet_sum, 2:3, function(x){
                const <- ifelse(max(x) == 0, 1, max(x))
                x / const
        } )


        # sp <- apply( beta.array.unet_sum, 1, max)

        # END sp ------------------------------------------------------------------

        params <- list(
                x = x,
                y = y,
                family = family,
                seq.alpha = seq.alpha,
                seq.df = seq.df,
                K = K,
                psub = psub,
                penalty.factor = penalty.factor,
                standardize.response = standardize.response,
                ...
        )


        return(list(sp = beta.array.unet_sum, sp01 = beta.array.unet_df, MaxConst = MaxConst, params = params))

}
