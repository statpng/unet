#' @name threshold
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
#'
#' Ycont <- Y[,c(1,4)]
#' Family_cont <-c("gaussian","gaussian")
#' penalty.factor <- c(0, rep(1, p))
#'
#'
#' lambda.list.mnet <- grid.lambda(x = X, y = Ycont,
#'                                 family = Family_cont,
#'                                 method = "mnet",
#'                                 iter = 10,
#'                                 seq.alpha = 1:9*0.1, n.lambda = 10,
#'                                 penalty.factor=penalty.factor)
#' lambda.vec.mnet <- lapply( lambda.list.mnet, function(x) seq(median(x), max(x), length.out=10) )
#'
#'
#' set.seed(1)
#' sp.mnet <- mnet(x = cbind(Z,X), y = Ycont,
#'                 seq.alpha = alpha.vec,
#'                 seq.lambda = lambda.vec.mnet, K = 10,
#'                 penalty.factor = penalty.factor )
#'
#'
#' set.seed(1)
#' sp.total_cont <- main( x = cbind(Z,X), y = Ycont,
#'                        family = Family_cont ,
#'                        method = c("enet", "mnet","unet","unet_df"),
#'                        seq.alpha = alpha.vec,
#'                        seq.lambda = list(enet=lambda.vec.enet[c(1,4)],
#'                                          mnet=lambda.vec.mnet),
#'                        seq.df = df.vec,
#'                        K = 10,
#'                        penalty.factor = penalty.factor )
#'
#'
#' head(sp.mnet$sp$mnet)
#' head(sp.total_cont$sp$unet)
#'
#' SelProb <- sp.total$sp
#' Threshold <- threshold(sp.total$params, nperm=10, seed=1)
#'
#' which( SelProb$enet[,1] > Threshold$empirical[,1,"enet"][10] )
#'
#'
#' @export threshold
threshold <- function(params, nperm = 100, seed = 1, const=NULL, perm.type=c("row-wise", "col-wise")) {

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

        n <- nrow( params$x )
        p <- sum( params$penalty.factor == 1 )
        q <- ncol( params$y )

        n.methods <- length(params$method)


        n.alpha <- length( params$seq.alpha )
        n.lambda <- length( params$seq.lambda )
        n.df <- length( params$seq.df )

        # empirical.result <- as.list(1:nperm)

        Empirical_Threshold <- array(0, dim=c(p, q, n.methods) )

        dimnames(Empirical_Threshold) <- list( 1:p, paste0("y.", 1:q), params$method )

        print("Permutation get started!")
        perm.params <- params
        pb <- txtProgressBar(min=0, max=nperm, style=3)
        for( perm.i in 1:nperm ){
                # set.seed((nperm+1)*j+perm.i)
                setTxtProgressBar(pb, perm.i)

                # permutation -------------------------------------------------------------
                set.seed( seed + perm.i )
                if( type == "row-wise" ){
                        perm.params$y <- params$y[sample(n),]
                } else if( type == "col-wise" ) {
                        perm.params$y <- apply( params$y, 2, function(yk) yk[sample(n)] )
                }

                fit.perm <- do.call("main", perm.params)

                # empirical.result[[perm.i]] <- fit.perm$sp

                # Empirical threshold -----------------------------------------------------
                for( mm in 1:n.methods ){
                        for( fdfd in 1:p ){
                                Empirical_Threshold[fdfd, , mm] <-
                                        Empirical_Threshold[fdfd, , mm] +
                                        apply( fit.perm$sp[[mm]], 2, function(x) sort(x, decreasing=TRUE)[ fdfd ] ) / nperm
                        }
                }
        }

        out <- as.list(1:n.methods)
        names(out) <- dimnames(Empirical_Threshold)[[3]]

        for( k in 1:n.methods ){
                out[[k]] <- Empirical_Threshold[,,k]
        }


        return( out )

}









threshold_local <- function(params, nperm = 100, seed = 1, const=NULL) {

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

        n <- nrow( params$x )
        p <- sum( params$penalty.factor == 1 )
        q <- ncol( params$y )

        n.methods <- length(params$method)


        n.alpha <- length( params$seq.alpha )
        n.lambda <- length( params$seq.lambda )
        n.df <- length( params$seq.df )

        # empirical.result <- as.list(1:nperm)

        Empirical_Threshold <- array(0, dim=c(p, q, n.methods) )

        dimnames(Empirical_Threshold) <- list( 1:p, paste0("y.", 1:q), params$method )

        print("Permutation get started!")
        perm.params <- params
        pb <- txtProgressBar(min=0, max=nperm, style=3)
        for( perm.i in 1:nperm ){
                # set.seed((nperm+1)*j+perm.i)
                setTxtProgressBar(pb, perm.i)

                # permutation -------------------------------------------------------------
                set.seed( seed + perm.i )
                perm.params$y <- params$y[sample(n),]

                fit.perm <- do.call("main", perm.params)

                # empirical.result[[perm.i]] <- fit.perm$sp

                # Empirical threshold -----------------------------------------------------
                for( mm in 1:n.methods ){
                        for( fdfd in 1:p ){
                                Empirical_Threshold[fdfd, , mm] <-
                                        Empirical_Threshold[fdfd, , mm] +
                                        apply( fit.perm$sp[[mm]], 2, function(x) sort(x, decreasing=TRUE)[ fdfd ] ) / nperm
                        }
                }
        }

        out <- as.list(1:n.methods)
        names(out) <- dimnames(Empirical_Threshold)[[3]]

        for( k in 1:n.methods ){
                out[[k]] <- Empirical_Threshold[,,k]
        }


        return( out )

}
