#### Sampler one non-binary matrix


"IPA.sampler.One_Non_Binary_Matrix" <- function(P, Iso, Int, no.its = 1100, burn = 100, delta = 1, allsamp = F, ratio.toll = 0.8, allsampcomp = NULL, v = TRUE, IT = 500) {
    
    
    Nc <- nrow(Iso)
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
            
            ### counting isotopes
            tmp <- matrix(Iso[sampcomp[-thism], ], ncol = Nc) * (Int[thism]/Int[-thism])
            ind.ones <- which((tmp >= ratio.toll) & (tmp <= (1/ratio.toll)))
            tmp[ind.ones] <- 1
            tmp[tmp != 1] <- 0
            p.iso <- colSums(tmp)
            
            
            ## normalising with deltas
            
            p.iso <- (p.iso + delta)/sum(p.iso + delta)
            
            ## merging scores
            po <- p.iso * P[thism, ]
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
