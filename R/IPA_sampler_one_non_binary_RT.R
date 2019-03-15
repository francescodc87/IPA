#### Sampler one non-binary matrix - Filter RT


"IPA.sampler.One_Non_Binary_Matrix_RT" <- function(P, Iso, Int, RT, RT.win = 3, ratio.toll = 0.8, no.its = 1100, burn = 100, delta = 1, allsamp = F, allsampcomp = NULL, 
    v = TRUE, IT = 500) {
    
    counter <- 0
    ind.rem <- lapply(RT, function(x) {
        counter <<- counter + 1
        checking.RT(RT, counter, RT.win)
    })
    
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
            tmp <- matrix(Iso[sampcomp[-ind.rem[[thism]]], ], ncol = Nc) * (Int[thism]/Int[-ind.rem[[thism]]])
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
