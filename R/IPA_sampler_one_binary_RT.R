#### Sampler one binary matrix - Filtering based on RT - to be used only for ADD and ISO binary

"IPA.sampler.One_Binary_Matrix_RT" <- function(P, Add, RT, RT.win = 3, no.its = 1100, burn = 100, delta = 1, allsamp = F, allsampcomp = NULL, v = TRUE, IT = 500) {
    
    
    counter <- 0
    ind.rem <- lapply(RT, function(x) {
        counter <<- counter + 1
        checking.RT(RT, counter, RT.win)
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
            
            ## normalising with deltas
            p.add <- (p.add + delta)/sum(p.add + delta)
            
            ## merging scores
            po <- p.add * P[thism, ]
            po <- po/sum(po)
            sampcomp[thism] <- multsample(po)
        }
        allsampcomp[it + olds, ] <- sampcomp
        if (v) {
            if (it%%100 == 0) {
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
