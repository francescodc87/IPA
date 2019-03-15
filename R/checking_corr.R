"checking.corr" <- function(Corr.vec, i, corr.thr) {
    out <- which(Corr.vec < corr.thr)
    out
}
