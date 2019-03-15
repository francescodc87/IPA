"IPA.sampler.Iso.Add" <- function(P, Iso, Add, Int, no.its = 1100, burn = 100, delta.iso = 0.5, delta.add = 0.5, allsamp = F, ratio.toll = 0.8, allsampcomp = NULL, v = TRUE, 
    IT = 500) {
    
    
    
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
    
    pot.add <- apply(Add[sampcomp, ], 2, sum)
    
    for (it in 1:no.its) {
        ordine <- sample(M)  # randomising the order used to cheack all assignments
        for (thism in ordine) {
            
            # counting adducts
            p.add <- pot.add - colSums(matrix(Add[sampcomp[thism], ], ncol = Nc))
            ### counting isotopes
            tmp <- matrix(Iso[sampcomp[-thism], ], ncol = Nc) * (Int[thism]/Int[-thism])
            ind.ones <- which((tmp >= ratio.toll) & (tmp <= (1/ratio.toll)))
            tmp[ind.ones] <- 1
            tmp[tmp != 1] <- 0
            p.iso <- colSums(tmp)
            
            
            ## normalising with deltas
            p.add <- (p.add + delta.add)/sum(p.add + delta.add)
            p.iso <- (p.iso + delta.iso)/sum(p.iso + delta.iso)
            
            ## merging scores
            po <- p.add * p.iso * P[thism, ]
            po <- po/sum(po)
            oldval <- sampcomp[thism]
            sampcomp[thism] <- multsample(po)
            if (oldval != sampcomp[thism]) {
                pot.add <- pot.add - Add[, oldval] + Add[, sampcomp[thism]]
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
