#' @import glmnet
#' @export grid.lambda
grid.lambda <-
  function(x, y,
           family,
           method,
           iter = 10,
           seq.alpha = NULL,
           n.lambda = NULL,
           psub = 0.5,
           penalty.factor = NULL,
           ...) {

    if( NROW(y) != nrow(x) ) stop("x and y should be equal length of row")

    if (is.null(seq.alpha)) seq.alpha <- 1:9 * 0.1
    if (is.null(n.lambda)) n.lambda <- 10
    if( is.null(penalty.factor)) penalty.factor <- rep(1, ncol(x))


    x <- as.matrix(x)
    n <- nrow(x)
    p <- ncol(x)
    nsub <- n * psub

    y.column.set <- seq_len(ncol(y))

    if( method %in% c("enet","unet") ) {
      lambda.array <- array( NA, dim = c(iter, n.lambda, length(seq.alpha), ncol(y), 2),
                             dimnames = list(
                               paste0("iter=", 1:iter),
                               paste0("lambda", 1:n.lambda),
                               paste0("alpha=", seq.alpha),
                               paste0("Y", 1:ncol(y)),
                               c("lambda", "df")
                             )
      )

      seq.lambda <- as.list(1:length(y.column.set))

      for (colcol in 1:length(y.column.set)) {
        FAMILY <- family[colcol]
        h <- unlist(y.column.set[colcol])

        lambda.vec <- NULL
        for (i in 1:iter) {
          for (j in 1:length(seq.alpha)) {
            wsub <- sample(n, nsub)
            xsub <- x[wsub, , drop = F]
            ysub <- y[wsub, , drop = F]
            if (FAMILY == "multinomial") {
              fitsub <- glmnet(
                x = xsub,
                y = y.multinom(ysub[,h,drop=T]),
                alpha = seq.alpha[j],
                family = FAMILY,
                nlambda = n.lambda,
                type.multinomial = "grouped",
                ... )
            } else {
              fitsub <- glmnet(
                x = xsub,
                y = ysub[, h],
                alpha = seq.alpha[j],
                family = FAMILY,
                nlambda = n.lambda,
                ...
              )
            }


            lambda.array[i, 1:length(fitsub$df), j, colcol, 1] <- fitsub$lambda
            lambda.array[i, 1:length(fitsub$df), j, colcol, 2] <- fitsub$df

            lambda.vec <- c(lambda.vec, fitsub$lambda)
          }
        }

        seq.lambda[[colcol]] <- lambda.vec

      }

    } else if (method == "mnet") {
      seq.lambda <- as.list(1:length(y.column.set))

      lambda.vec <- NULL
      for (i in 1:iter) {
        for (j in 1:length(seq.alpha)) {
          wsub <- sample(n, nsub)
          xsub <- x[wsub, , drop = F]
          ysub <- y[wsub, , drop = F]
          fitsub <-
            glmnet(
              x = xsub,
              y = ysub,
              alpha = seq.alpha[j],
              family = "mgaussian",
              nlambda = n.lambda,
              ...
            )

          lambda.vec <- c(lambda.vec, fitsub$lambda)
        }
      }

      for (colcol in 1:length(y.column.set)) {
        seq.lambda[[colcol]] <- lambda.vec

      }

    }

    return(seq.lambda)

  }
