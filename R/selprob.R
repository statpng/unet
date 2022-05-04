selprob <- function(array, K = 100) {
  Margin.rep <- which(!names(dimnames(array)) %in% c("Replications"))
  count.array <- apply(array, Margin.rep, mean)
  Margin.tuning <- which(!names(dimnames(count.array)) %in% c("Alpha", "Lambda"))
  out.sp <- as.matrix(apply(count.array, Margin.tuning, max))

  c("Variable", "Alpha", "Lambda", "Phenotype", "Replications")

  Margin.phenotype <- which(names(dimnames(count.array)) %in% c("Phenotype"))

  PI.seq <- c(1, 5, 10, 50, 100)
  threshold <- matrix(NA, length(PI.seq), dim(count.array)[Margin.phenotype])
  dimnames(threshold) <- list(PI.seq, paste0("y", seq_len(ncol(threshold))))
  for (colcol in 1:ncol(threshold)) {
    qhat <- (sum(count.array[, , , colcol]) * K) / (prod(dim(count.array)[-Margin.tuning]) * K)
    for (pipi in 1:length(PI.seq)) {
      PI <- PI.seq[pipi]
      threshold[pipi, colcol] <- min(qhat ^ 2 / (2 * PI * dim(count.array)[1]) + 1 / 2, 1)
    }
  }

  list(sp = out.sp, threshold = threshold)
}
