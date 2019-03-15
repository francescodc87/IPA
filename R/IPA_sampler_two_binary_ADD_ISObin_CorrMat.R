### two binary matrices, filtering RT for ADD and ISO BIO, Bio


"IPA.sampler.Two_Binary_Matrix_CorrMat_Add_IsoBin" <- function(P, Add, Iso, Corr.matrix, corr.thr = 0.8, no.its = 1100, burn = 100, delta.add = 0.5, delta.iso = 0.1, allsamp = F, 
    allsampcomp = NULL, v = TRUE, IT = 500) {
    
    
    counter <- 0
    ind.rem <- lapply(rep(NA, nrow(P)), function(x) {
        counter <<- counter + 1
        checking.corr(Corr.matrix[, counter], counter, corr.thr)
    })
    
    Nc <- nrow(Add)
    M <- nrow(P)
    
    if (is.null(allsampcomp)) {
        sampcomp <- apply(P, 1, multsample)
        allsampcomp <- matrix(0, no.its, M)
        olds <- 0
    } else {
        olds <- nrow(allsampcomp)
        sampcomp <- allsampcomp[olds, ]
        allsampcomp <- rbind(allsampcomp, matrix(0, no.its, M))
    }
    
    
    for (it in 1:no.its) {
        ordine <- sample(M)  # randomising the order used to cheack all assignments
        for (thism in ordine) {
            
            # counting adducts
            p.add <- colSums(matrix(Add[sampcomp[-ind.rem[[thism]]], ], ncol = Nc))
            ### counting isotopes
            p.iso <- colSums(matrix(Iso[sampcomp[-ind.rem[[thism]]], ], ncol = Nc))
            
            
            ## normalising with deltas
            p.add <- (p.add + delta.add)/sum(p.add + delta.add)
            p.iso <- (p.iso + delta.iso)/sum(p.iso + delta.iso)
            
            ## merging scores
            po <- p.add * p.iso * P[thism, ]
            po <- po/sum(po)
            sampcomp[thism] <- multsample(po)
        }
        allsampcomp[it + olds, ] <- sampcomp
        if (v) {
            if (it%%IT == 0) {
                # Print on the screen some message
                cat(paste0(round((it * 100)/no.its, 1), "%", "\n"))
            }
        }
    }
    
    
    post <- t(apply(allsampcomp, 2, compute.post, burn = burn, no.its = no.its + olds, Nc = Nc))
    if (allsamp) {
        out <- list(Post = post, allsampcomp = allsampcomp)
        return(out)
    } else {
        return(post)
    }
}
