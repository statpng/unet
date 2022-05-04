# x <- rbind( mnormt::rmnorm(100, rep(0, 500), diag(500)), mnormt::rmnorm(100, rep(1, 500), diag(500)) )
# y <- data.frame(y1=rep(0:1, each=100), y2=x%*%rep(0.1,500)+rnorm(200), y3=rbinom(200, 5, prob=rep(c(0.4,0.6), each=100)))
# family <- rep("gaussian", ncol(y))
# type.multivariate <- "unet"
# seq.alpha <- 1*0.1
# seq.lambda <- replicate(5, seq(0.5, 6.0, length.out=10), simplify = F)
# K <- 100
# psub <- 0.5


y.multinom <- function(yk){
  model.matrix(~ . - 1, data=as.data.frame(as.factor(yk) %>% droplevels()))
}


th <- function(q, alpha, p)
  q ^ 2 / (2 * alpha * p) + 1 / 2

varcov_rho <- function(p, rho) {
  out <- matrix(rho, p, p)
  diag(out) <- 1
  return(out)
}

png.varcov <-
  function(p,
           rho = 0,
           type = NULL,
           Omega = 0.001,
           PropOfNeg = 0.25) {
    # Depends : varcov_rho

    if (!type %in% c("arcov", "random_sign", "all.equal"))
      stop("type should be one among arcov, random_sign, and all.equal")

    if (is.null(type))
      stop("'type' should be entered")

    if (type == "arcov") {
      out <- outer(1:p, 1:p, function(x, y)
        rho ^ abs(x - y))
      return(out)
    }

    if (type == "random_sign") {
      if (PropOfNeg < 0 |
          PropOfNeg > 0.5)
        stop("PropOfNeg must be in [0,0.5].")


      e.rho <-
        replicate(p * (p - 1) / 2, runif(1, 0, 1) * (-1) ^ rbinom(1, 1, PropOfNeg))
      e.varcov <- matrix(1, p, p)
      e.varcov[upper.tri(e.varcov)] <- e.rho
      e.varcov[lower.tri(e.varcov)] <-
        t(e.varcov)[lower.tri(e.varcov)]

      if (isSymmetric(e.varcov)) {
        x <- (e.varcov %*% e.varcov) + diag(Omega, p)
        e.varcov2 <- (x / sqrt(diag(x) %*% t(diag(x))))
        return(e.varcov2)
      } else {
        stop("Error")
      }
    }

    if (type == "all.equal") {
      # if( rho < 0 ) stop("rho have to be equal or greater than 0")
      e.varcov2 <- varcov_rho(p = p, rho = rho)
      return(e.varcov2)
    }

    # e.varcov2 %>% .[upper.tri(.)] %>% hist(main="Error-term correlation", xlab=expression(rho,"e"))+
  }

png.msnp <-
  function(n,
           p,
           M = 25,
           snp.rho = 0.6,
           y.rho = 0.3,
           maf.min = 0.01,
           sp.mu,
           cp.mu,
           ncp,
           nsp,
           ncp.size,
           cov.type = "all.equal",
           y.p.neg = NULL,
           omega = NULL,
           standardization = TRUE) {
    # Depends: png.snp, png.varcov, mnormt

    # sp.mu <- rnorm(1, 0.3, 0.05)
    # cp.mu <- rnorm(1, 0.3, 0.05)

    PNG.SNP <- png.snp(
      n = n,
      p = p,
      rho = snp.rho,
      min.maf = maf.min
    )
    if (is.matrix(y.rho)) {
      VAR <- y.rho
    } else {
      VAR <- png.varcov(M, rho = y.rho, type = cov.type)
    }
    SNP <- PNG.SNP$snp
    dimnames(SNP) <- list(paste0("N", 1:n), paste0("snp", 1:p))
    MAF <- PNG.SNP$MAF

    nsig <- ncp + nsp
    # SNP <- do.call( "cbind", lapply(MAF, function(x) rbinom(n, 2, x) ) )
    ##  range( cor(SNP)[upper.tri(matrix( 0, ncol(SNP), ncol(SNP) ) )] )

    if (nsp %% M != 0)
      warnings("The number of single-phenotypic variants is not proportional to M!")
    if (ncp > 0)
      var.cp <- (1:nsig)[(1:nsig) %in% seq_len(ncp)]
    if (nsp > 0)
      var.sp <- (1:nsig)[!(1:nsig) %in% seq_len(ncp)]

    beta <- replicate(M, rep(0, p))
    if (ncp > 0) {
      for (m in 1:ncp.size) {
        beta[var.cp * 10, m] <- cp.mu
      }
    }
    if (nsp > 0)  {
      var.sp.list <-
        tapply(var.sp, rep(1:M, each = ceiling(nsp / M))[1:nsp], list)
      for (m in 1:length(var.sp.list)) {
        beta[var.sp.list[[m]] * 10, m] <- sp.mu
      }
    }


    #  for( ss in 1:nsig ) {
    #	if( ss %in% var.cp ){
    #		for( m in 1:ncp.size ){
    #			beta[ var.cp[ss]*10, m] <- cp.mu
    #		}
    #	} else if ( ss %in% var.sp ){
    #		for( m in 1:length(var.sp) ){
    #			beta[ var.sp[m]*10, ((m-1)%%M+1)] <- sp.mu
    #		}
    #	}
    # }


    true <- which(beta != 0, arr.ind = TRUE)
    if (ncp > 0)  {
      true.cp <- unique(true[duplicated(true[, 1]), 1, drop = T])
    } else {
      true.cp <- NULL
    }
    if (nsp > 0)  {
      if (is.null(true.cp)) {
        true.sp <- true[, 1]
      } else {
        true.sp <- true[!true[, 1] %in% true.cp, 1]
      }
      true.sp <-
        tapply(true.sp , rep(1:M, each = ceiling(nsp / M))[1:nsp], list)
    } else {
      true.sp <- NULL
    }



    Y <- SNP %*% beta + rmnorm(n = n, varcov = VAR)
    #  cor(Y)
    #  corrplot::corrplot(cor(y), tl.pos = "n")

    ValuesOfArguments =
      list(
        n = n,
        p = p,
        M = M,
        snp.rho = snp.rho,
        y.rho = y.rho,
        maf.min = maf.min,
        sp.mu = sp.mu,
        cp.mu = cp.mu,
        ncp = ncp,
        nsp = nsp,
        cov.type = cov.type
      )

    if (standardization)
      SNP <- scale(SNP)

    Data <-
      list(
        snp = SNP,
        y = Y,
        maf = MAF,
        beta = beta,
        true = true,
        true.cp = true.cp,
        true.sp = true.sp,
        args = ValuesOfArguments
      )
    return(Data)
  }


png.snp <-
  function(n,
           p,
           rho = 0.95,
           min.maf = 0.05,
           maf.type = c("unif", "beta", "chisq")) {
    X <-
      do.call("cbind", lapply(1:20, function(x)
        mnormt::rmnorm(n, varcov = png.varcov(
          p = (p / 20),
          rho = rho,
          type = "arcov"
        ))))
    Y <-
      do.call("cbind", lapply(1:20, function(x)
        mnormt::rmnorm(n, varcov = png.varcov(
          p = (p / 20),
          rho = rho,
          type = "arcov"
        ))))

    if (maf.type == "unif") {
      MAF <- runif(p, min.maf, 0.5)
      # ref
      # Dai, M., Ming, J., Cai, M., Liu, J., Yang, C., Wan, X., & Xu, Z. (2017). IGESS: a statistical approach to integrating individual-level genotype data and summary statistics in genome-wide association studies. Bioinformatics, 33(18), 2882-2889.
      # Yang, Y., Shi, X., Jiao, Y., Huang, J., Chen, M., Zhou, X., ... & Liu, J. (2020). CoMM-S2: a collaborative mixed model using summary statistics in transcriptome-wide association studies. Bioinformatics, 36(7), 2009-2016.
      # Jiang, W., & Yu, W. (2017). Controlling the joint local false discovery rate is more powerful than meta-analysis methods in joint analysis of summary statistics from multiple genome-wide association studies. Bioinformatics, 33(4), 500-507.
    } else if (maf.type == "beta") {
      MAF <-
        rbeta(p, 0.14, 0.73) %>% {
          ifelse(. < 0.5, ., 1 - .) * (1 - min.maf * 2) + min.maf
        }
      # ref
      # Ionita-Laza, I., Lange, C., & Laird, N. M. (2009). Estimating the number of unseen variants in the human genome. Proceedings of the National Academy of Sciences, 106(13), 5008-5013.
    } else if (maf.type == "chisq") {
      MAF <- truncnorm::rtruncnorm(
        p,
        a = -0.65,
        b = 0.65,
        mean = 0,
        sd = 5
      ) ^ 2
      # MAF <- truncnorm::rtruncnorm(p, a=sqrt(min.maf), b=0.65, mean=0, sd=5)^2
      # No reference
    }


    for (j in seq_len(p)) {
      perct_X <- rank((X[, j])) / length(X[, j])
      perct_Y <- rank((Y[, j])) / length(Y[, j])

      X[perct_X <= MAF[j], j] <- 1
      X[perct_X >  MAF[j], j] <- 0
      Y[perct_Y <= MAF[j], j] <- 1
      Y[perct_Y >  MAF[j], j] <- 0
    }

    Data <- rbind(X + Y)

    return(list(snp = Data, MAF = MAF))
  }
