"compute.post" <- function(allsampcomp, no.its, burn, Nc) {
    tmp <- table(allsampcomp[(burn + 1):no.its])/(no.its - burn)
    out <- rep(0, Nc)
    out[as.numeric(names(tmp))] <- tmp
    out
}
