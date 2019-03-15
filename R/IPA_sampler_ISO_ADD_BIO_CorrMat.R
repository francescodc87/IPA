## ADD ISO BIO -- RT filtering

"IPA.sampler.Add.Iso.Bio.CorrMat" <- function(P, Add, Iso, Bio, Int, Corr.matrix, corr.thr = 0.8, no.its = 1100, burn = 100, delta.add = 0.5, delta.iso = 0.5, delta.bio = 1, 
    allsamp = F, ratio.toll = 0.8, allsampcomp = NULL, v = TRUE, IT = 500) {
    
    
    
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
    
    pot.bio <- apply(Bio[sampcomp, ], 2, sum)
    
    
    for (it in 1:no.its) {
        ordine <- sample(M)  # randomising the order used to cheack all assignments
        for (thism in ordine) {
            
            # counting adducts
            p.add <- colSums(matrix(Add[sampcomp[-ind.rem[[thism]]], ], ncol = Nc))
            ### counting isotopes
            tmp <- matrix(Iso[sampcomp[-ind.rem[[thism]]], ], ncol = Nc) * (Int[thism]/Int[-ind.rem[[thism]]])
            ind.ones <- which((tmp >= ratio.toll) & (tmp <= (1/ratio.toll)))
            tmp[ind.ones] <- 1
            tmp[tmp != 1] <- 0
            p.iso <- colSums(tmp)
            
            ## counting biotransformations
            p.bio <- pot.bio - colSums(matrix(Bio[sampcomp[thism], ], ncol = Nc))
            
            ## normalising with deltas
            p.add <- (p.add + delta.add)/sum(p.add + delta.add)
            p.iso <- (p.iso + delta.iso)/sum(p.iso + delta.iso)
            p.bio <- (p.bio + delta.bio)/sum(p.bio + delta.bio)
            
            ## merging scores
            po <- p.add * p.iso * p.bio * P[thism, ]
            po <- po/sum(po)
            oldval <- sampcomp[thism]
            sampcomp[thism] <- multsample(po)
            if (oldval != sampcomp[thism]) {
                pot.bio <- pot.bio - Bio[, oldval] + Bio[, sampcomp[thism]]
            }
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
