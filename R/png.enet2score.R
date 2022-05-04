#' @export png.enet2score
png.enet2score <- function(Array, delta = 1.0, weight = TRUE){
        p <- nrow(Array)
        q <- ncol(Array)


        Array.sort <- apply(Array, 2, sort, decreasing=TRUE)
        Array.order <- apply(Array, 2, order, decreasing=TRUE)
        Array.sum <- apply(Array, 2, sum)
        # Array.sort %>% head(10)
        # Array.order %>% head(10)
        # Array.sum %>% head(10)
        sum.min <- floor( min(Array.sum) * delta )
        sum.which.min <- which.min(Array.sum)

        dup.label <- function(vec){
                count <- 0
                out <- NULL
                for( lev in unique(vec) ){
                        count <- count + 1
                        out[which( vec == lev )] <- count
                }
        }

        Array.sort.label <- apply( Array.sort, 2, dup.label )

        wh.cumsum.zero <- minvalue.cumsum.zero <- NULL
        for( k in 1:q ){
                wh.cumsum.zero[k] <- min( which( cumsum( Array.sort[,k] ) >= sum.min ) )
                minvalue.cumsum.zero[k] <- cumsum( Array.sort[,k] )[ wh.cumsum.zero[k] ]

                # Array.sort[which( cumsum( Array.sort[,k] ) >= sum.min ), k]
                # cumsum( Array.sort[,k] )[which( cumsum( Array.sort[,k] ) >= sum.min )] %>% head

                if( sum( Array.sort.label[,k] == Array.sort.label[wh.cumsum.zero[k], k] ) > 1 ){

                        wh.cumsum.zero[k] <- max( which( Array.sort.label[,k] == Array.sort.label[wh.cumsum.zero[k], k] ) )
                        minvalue.cumsum.zero[k] <- cumsum( Array.sort[,k] )[ wh.cumsum.zero[k] ]

                }

        }

        # wh.cumsum.zero
        # minvalue.cumsum.zero

        # Thresholding
        out <- NULL
        for( k in 1:q ){
                out[[k]] <- Array[ , k ][ Array.order[ 1:wh.cumsum.zero[k], k ], drop=F ]
        }


        out.sum <- matrix( 0, nrow = p, ncol = length(out) )
        for( k in 1:length(out) ){
                if( weight ){
                        out.sum[ as.numeric( names(out[[k]]) ), k ] <- out[[k]] * sum.min / max(1, minvalue.cumsum.zero[k])
                } else {
                        out.sum[ as.numeric( names(out[[k]]) ), k ] <- out[[k]]
                }

        }


        # -- Validation --
        # sapply( 1:ncol(Array), function(k){
        #   Array[ 1:20, k ] * sum.min / minvalue.cumsum.zero[k]
        # }) %>% apply(1, sum) %>% {./100}
        #
        # apply( out.sum, 1, sum )[1:20] / 100
        # -- Validation --



        return(
                list( sumsv = Array.sum,
                      sumsv.dlt = minvalue.cumsum.zero,
                      nvars.dlt = wh.cumsum.zero,
                      score = out.sum )
        )




}
