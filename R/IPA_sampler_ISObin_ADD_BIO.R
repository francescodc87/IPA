"IPA.sampler.Add.IsoBin.Bio" <- function(P, Add, Iso, Bio, no.its = 1100, burn = 100, delta.add = 0.5, delta.iso = 0.5, delta.bio = 1, allsamp = F, allsampcomp = NULL, v = TRUE, 
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
    
    pot.bio <- apply(Bio[sampcomp, ], 2, sum)
    pot.add <- apply(Add[sampcomp, ], 2, sum)
    pot.iso <- apply(Iso[sampcomp, ], 2, sum)
    
    
    for (it in 1:no.its) {
        ordine <- sample(M)  # randomising the order used to cheack all assignments
        for (thism in ordine) {
            
            # counting adducts
            p.add <- pot.add - colSums(matrix(Add[sampcomp[thism], ], ncol = Nc))
            ### counting isotopes
            p.iso <- pot.iso - colSums(matrix(Iso[sampcomp[thism], ], ncol = Nc))
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
                pot.add <- pot.add - Add[, oldval] + Add[, sampcomp[thism]]
                pot.iso <- pot.iso - Iso[, oldval] + Iso[, sampcomp[thism]]
                
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

