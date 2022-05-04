#' @name mnet
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
#' @return A data.frame with each column representing the selection probability for each phenotype
#' @import mnormt
#' @import glmnet
#' @export mnet
mnet <-
  function(x,
           y,
           seq.alpha = NULL,
           seq.lambda = NULL,
           K = 100,
           psub = 0.5,
           standardize.response = TRUE,
           penalty.factor = NULL,
           ...) {



    if( is.list(seq.lambda) & length(seq.lambda) != ncol(y) ) stop("The length of seq.lambda should be equal to ncol(y).")

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
    n.lambda <- unique( sapply(seq.lambda, length) )

    # IF sp -------------------------------------------------------------------
    beta.array.mnet <- array(0, dim = c(p, n.alpha, n.lambda, ncol(y)),
                             dimnames = list( paste0("", 1:p),  # paste0("v", 1:p),
                                              paste0("", seq.alpha),   # paste0("Alpha=",seq.alpha),
                                              paste0("", seq_len(n.lambda)), # paste0("Lambda=", seq_len(n.lambda)),
                                              paste0("", seq_len(ncol(y))) ) )

    names(attributes(beta.array.mnet)$dimnames) <- c("Variable", "Alpha", "Lambda", "Phenotype")


    for( kk in 1:K ){
      if( kk %% 10 == 0 ) cat("Iteration :", kk, "/", K, "\n")
      start <- proc.time();
      wsub <- sample(n, nsub)
      xsub <- x[wsub,,drop=F];
      ysub <- y[wsub,,drop=F];

      for( aa in 1:length(seq.alpha) ){

        mgaussian.fit <- glmnet(x=xsub,
                                y=ysub,
                                alpha=seq.alpha[aa],
                                lambda=unlist(seq.lambda[1]),
                                family="mgaussian",
                                standardize.response=standardize.response,
                                penalty.factor = penalty.factor, ...)$beta[[1]][wh.var,]
        for( colcol in 1:ncol(y) ){
            beta.array.mnet[,aa,,colcol] <- beta.array.mnet[,aa,,colcol] + as.numeric( mgaussian.fit != 0 )/K
        }

      }
      # end alpha ---------------------------------------------------------------

      end <- proc.time();
      # cat("In this iteration, the elaspsed time=", (end-start)[3], "\n");
      if(kk == 1) cat("The expected remained time is", (end - start)[3] * (K - kk), "\n")

    }

    # end REP -----------------------------------------------------------------

    # out <- list(mnet = slam::as.simple_sparse_array( beta.array.mnet ))
    out <- list(mnet = apply( beta.array.mnet, c(1,4), max ) )

    # END sp ------------------------------------------------------------------

    return(out)

  }
