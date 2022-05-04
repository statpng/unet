#' @name df2lambda
#' @title Selection probabilities of the regularization model using the elastic-net penalty with multiresponse gaussian family (MNET).
#' @description
#' The penalized regression model using elastic-net penalty with multiresponse gaussian family (MNET) aims to solve the following problem:
#'
#' @export df2lambda
df2lambda <- function(x, y,
                      alpha,
                      family, type.multinomial, seq.df,
                      penalty.factor = NULL, verbose=TRUE){

  if( is.null( penalty.factor ) ) penalty.factor <- rep(1, ncol(x))
  if( sum( penalty.factor==0 ) > min(seq.df) ) seq.df[which.min(seq.df)] <- seq.df[which.min(seq.df)] + sum( penalty.factor==0 )

  get_range_glmnet <- function(glmnet_fit, DFDF){

    df.loss <- abs(glmnet_fit$df - DFDF)
    out.df <- glmnet_fit$df[ which.min(df.loss) ]
    out.lambda <- glmnet_fit$lambda[ which.min(df.loss) ]

    if( sum( df.loss == 0 ) > 0 ){
      lambda_upper <- glmnet_fit$lambda[ median( which( df.loss == 0 ) ) ]
      lambda_lower <- glmnet_fit$lambda[ median( which( df.loss == 0 ) ) ]
    } else {
      lambda_minmax <- glmnet_fit$lambda[ order(df.loss, decreasing=FALSE)[1:2] ]
      lambda_upper <- max( lambda_minmax )
      lambda_lower <- min( lambda_minmax )
    }

    #
    # df_min <- glmnet_fit$df[which.min(df.loss)]
    # df_upper <- ifelse( all(!glmnet_fit$df > df_min),
    #                     max(glmnet_fit$df),
    #                     ifelse( which.min(df.loss)==length(df.loss),
    #                             length(df.loss),
    #                             min( glmnet_fit$df[ glmnet_fit$df > df_min ] ) ) )
    # df_lower <- ifelse( which.min(df.loss)==1, min(glmnet_fit$df), max( glmnet_fit$df[ glmnet_fit$df < df_min ] ) )
    #
    # lambda_min <- glmnet_fit$lambda[which.min(df.loss)]
    # lambda_lower <- min( glmnet_fit$lambda[glmnet_fit$df == df_upper] )
    # lambda_upper <- max( glmnet_fit$lambda[glmnet_fit$df == df_lower] )

    list( df = out.df, lambda = out.lambda, minmax=c(lower=lambda_lower, upper=lambda_upper) )
  }



  fit.lambda_range <- glmnet( x=x, y=y,
                              alpha=alpha, nlambda=nrow(x)*10,
                              family=family, type.multinomial = type.multinomial, dfmax = 2*nrow(x),
                              penalty.factor = penalty.factor)

  lambda.sequence <- fit.lambda_range$lambda
  df.sequence <- fit.lambda_range$df


  out.lambda.seq <- NULL
  for( DFDF in seq.df ){
    fit.lambda_minmax <- get_range_glmnet(fit.lambda_range, DFDF)
    avg.lambda.minmax <- mean(fit.lambda_minmax$minmax)

    out.lambda.seq <- c(out.lambda.seq, avg.lambda.minmax)
  }

  tmp_glmnet_fit <-
    glmnet( x=x,
            y=y,
            alpha=alpha,
            lambda=out.lambda.seq,
            family=family,
            type.multinomial = type.multinomial,
            penalty.factor = penalty.factor)

  if( verbose ){
    print( tmp_glmnet_fit )

    cat( paste0("Target d.f.= ", seq.df, ", Estimated d.f.= ", tmp_glmnet_fit$df, collapse = "\n"), "\n" )

  }

  list( DF=seq.df, df=tmp_glmnet_fit$df, lambda=out.lambda.seq )

}





# df2lambda <- function(x, y,
#                       alpha, nlambda,
#                       family, type.multinomial, DFDF, tol, ...){
#
#   fit.lambda_range <- glmnet( x=x, y=y,
#                               alpha=alpha, nlambda=nlambda, dfmax=DFDF+1,
#                               family=family, type.multinomial = type.multinomial, ...)
#
#     df.loss <- abs(fit.lambda_range$df - DFDF)
#   df_min <- fit.lambda_range$df[which.min(df.loss)]
#   df_upper <- ifelse( all(!fit.lambda_range$df > df_min),
#                       max(fit.lambda_range$df),
#                       ifelse( which.min(df.loss)==length(df.loss),
#                               length(df.loss),
#                               min( fit.lambda_range$df[ fit.lambda_range$df > df_min ] ) ) )
#   df_lower <- ifelse( which.min(df.loss)==1, 0, max( fit.lambda_range$df[ fit.lambda_range$df < df_min ] ) )
#
#   lambda_min <- fit.lambda_range$lambda[which.min(df.loss)]
#   lambda_upper <- min( fit.lambda_range$lambda[fit.lambda_range$df == df_upper] )
#   lambda_lower <- max( fit.lambda_range$lambda[fit.lambda_range$df == df_lower] )
#
#   fit.opt <- optimize( function(LAMBDA, DF, ...){
#     fit2 <- glmnet( ..., lambda=LAMBDA );
#     abs(fit2$df - DF)},
#     interval = c(lambda_lower, lambda_upper),
#     DF=DFDF, x=x, y=y, alpha=alpha, family=family, type.multinomial=type.multinomial, dfmax=DFDF+1,
#     tol=tol )
#
#   opt.lambda=fit.opt$minimum
#   opt.lambda
#
# }





#' @Example
#' set.seed(1)
#' n=100; p=1000; q=4
#'
#' X <- replicate(p, rbinom(n,2,0.2)) #generate the SNP X=0,1,2
#' b <- matrix(0, p, q)
#' b[1:5,1:2] <- 1.0
#'
#' Z <- replicate(1, rnorm(n))
#' g <- matrix(0, 1, q)
#' g[1,1:4] <- 0.1
#'
#' Y <- cbind(Z, X)%*%rbind(g, b) + replicate(q, rnorm(n,1))
#' Y[,2:3] <- apply(Y[,2:3], 2, function(yk) ifelse(yk > median(yk), 1, 0) )
#' Family<-c("gaussian","binomial","binomial","gaussian")
#' penalty.factor <- c(0, rep(1, p))

# library(glmnet)
# df2lambda(x = cbind(Z, X), y = Y[,1], alpha=1, family=Family[1], seq.df=1:20*5, penalty.factor=penalty.factor, verbose=TRUE)
# glmnet(x = cbind(Z, X), y = Y[,1], alpha=1, lambda = 0.346986, family=Family[1], seq.df=10, penalty.factor=penalty.factor)



# x = cbind(Z, X)
# y= Y[,1]
# alpha=1
# family=Family[1]
# seq.df= 81:90
# penalty.factor=penalty.factor
# verbose=TRUE
#
#
#
# fit.lambda_range <- glmnet( x=x, y=y,
#                             alpha=alpha, nlambda=nrow(x)*50,
#                             family=family, type.multinomial = type.multinomial, dfmax = 2*nrow(x),
#                             penalty.factor = penalty.factor)
#
# lambda.sequence <- fit.lambda_range$lambda
# df.sequence <- fit.lambda_range$df
#
#
# out.lambda.seq <- NULL
# for( DFDF in seq.df ){
#   fit.lambda_minmax <- get_range_glmnet(fit.lambda_range, DFDF)
#   avg.lambda.minmax <- mean(fit.lambda_minmax$minmax)
#
#   out.lambda.seq <- c(out.lambda.seq, avg.lambda.minmax)
# }
#
# tmp_glmnet_fit <-
#   glmnet( x=x,
#           y=y,
#           alpha=alpha,
#           lambda=out.lambda.seq,
#           family=family,
#           type.multinomial = type.multinomial,
#           penalty.factor = penalty.factor)
#
# if( verbose ){
#   print( tmp_glmnet_fit )
#
#   cat( paste0("Target d.f.= ", seq.df, ", Estimated d.f.= ", tmp_glmnet_fit$df, collapse = "\n"), "\n" )
#
# }
#
# list( DF=seq.df, df=tmp_glmnet_fit$df, lambda=out.lambda.seq )
