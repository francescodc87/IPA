#### Sampler one binary matrix

"IPA.sampler.One_Binary_Matrix" <- function(P, Bio, no.its = 1100, burn = 100, delta = 1, allsamp = F, allsampcomp = NULL, v = TRUE, IT = 500) {
    
    Nc <- nrow(Bio)
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
            ## counting biotransformations
            p.bio <- pot.bio - colSums(matrix(Bio[sampcomp[thism], ], ncol = Nc))
            
            ## normalising with deltas
            p.bio <- (p.bio + delta)/sum(p.bio + delta)
            
            ## merging scores
            po <- p.bio * P[thism, ]
            po <- po/sum(po)
            oldval <- sampcomp[thism]
            sampcomp[thism] <- multsample(po)
            if (oldval != sampcomp[thism]) {
                pot.bio <- pot.bio - Bio[, oldval]
                pot.bio <- pot.bio + Bio[, sampcomp[thism]]
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
