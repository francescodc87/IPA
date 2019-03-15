"checking.rel.id" <- function(rel.id, i) {
    out <- c(i, which(rel.id != rel.id[i]))
    out
}
