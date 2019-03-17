#' @title Computing prior probabilities associated with putative annotations
#'
#' @description
#' This functions takes as input the output of the find.hits() function and computes the prior
#' probabilities for the putative annotations
#'
#'
#' @param Hits  The output of find.hits() function
#' @param dataset A matrix containing the measured data, organized in 3 colums: mz, RT and Int
#' @param pk A vector of length nrow(Hits$all.formulas). Values between 0 and 1 expressing initial confidence of the presence of each formula
#' @param ppm A number indicating the instrument accuracy to be considered
#' @param unknown.ppm The ppm number to be assigned to unknown (default NA - No unknown is considered)
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 120)
#' @param pr.lim A number inidicating the lowest probability value allowed (default 1e-05)
#'
#' @return A list containing the matrix of the prior probabilities and the id.masses and all.formulas objects
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso find.Hits IPAposteriors
#'
#' @import Matrix
#' @import enviPat
#' @import stringr
#' @importFrom utils combn
#' @export


### the funtion needs enviPat
"compute.Priors" <- function(Hits, dataset, pk = rep(1, nrow(Hits$all.formulas)), ppm, unknown.ppm = NA, v = T, IT = 120, pr.lim = 1e-05) {
    cat("Computing Priors... \n")
    # defing the number of masses and the number of the compounds
    compounds.mass <- as.numeric(Hits$all.formulas[, 5])
    Nc <- nrow(Hits$all.formulas)
    mass <- as.numeric(dataset[Hits$id.masses, 1])
    M <- length(mass)
    ### evaluating precision
    deltaMs <- ppm * mass * (1e-06)
    sigma <- deltaMs/2
    precision <- 1/(sigma^2)
    rm(sigma, deltaMs)
    # evaluation of prior probabilities (likelihood based only on mass) initialize some variables

    if (!is.na(unknown.ppm)) {
        pr <- Matrix(0, M, (Nc + 1))
        for (k in 1:M) {
            pr[k, 1:Nc] <- ((exp((-0.5 * precision[k]) * ((compounds.mass - mass[k])^2))) * pk)
            delta.unknown <- unknown.ppm * mass[k] * 1e-06
            pr[k, Nc + 1] <- ((exp((-0.5 * precision[k]) * ((delta.unknown)^2))))
            if (v) {
                if (k%%IT == 0) {
                  # Print on the screen some message
                  cat(paste0(round((k * 100)/M, 1), "%", "\n"))
                }
            }
        }
        all.formulas1 = rbind(Hits$all.formulas, c("unknown", rep(NA, 7)))
    } else {
        pr <- Matrix(0, M, Nc)
        for (k in 1:M) {
            pr[k, 1:Nc] <- ((exp((-0.5 * precision[k]) * ((compounds.mass - mass[k])^2))) * pk)
            if (v) {
                if (k%%IT == 0) {
                  # Print on the screen some message
                  cat(paste0(round((k * 100)/M, 1), "%", "\n"))
                }
            }
        }
        all.formulas1 = Hits$all.formulas
    }
    for (k in 1:nrow(pr)) {
        pr[k, ] <- pr[k, ]/sum(pr[k, ])
    }
    idx.lim <- which(as.matrix(pr)<pr.lim)
    pr[idx.lim] <- 0
    cat("\n", length(idx.lim))
    ind.C <- which(apply(pr, 2, sum, na.rm = T) > 0)
    all.formulas = all.formulas1[ind.C, ]  ### of course it is considering the last column even if there isn't the unknown entry in Hits!!
    if (length(ind.C) < ncol(pr)) {
        ## I have to check for orphans isotopes
        cat("Checking for orphan isotopes... \n")
        IDs <- unique(all.formulas[, 1:2])
        orphans <- NULL
        for (k in 1:nrow(IDs)) {
            ind <- which(all.formulas[, 1] == IDs[k, 1] & all.formulas[, 7] == "mono" & all.formulas[, 2] == IDs[k, 2])
            if (length(ind) == 0) {
                orphans <- c(orphans, which(all.formulas1[, 1] == IDs[k, 1] & all.formulas1[, 2] == IDs[k, 2]))
            }
            if (v) {
                if (k%%IT == 0) {
                  # Print on the screen some message
                  cat(paste0(round((k * 100)/nrow(IDs), 1), "%", "\n"))
                }
            }
        }
    }
    if (!is.null(orphans)) {
        ind.C <- ind.C[which(!ind.C %in% orphans)]
    }
    pr <- pr[, ind.C]
    ind.M <- which(apply(pr[,1:(nrow(pr)-1)], 1, sum, na.rm = T) > 0)
    pr <- pr[ind.M, ]
    all.formulas = all.formulas1[ind.C, ]
    id.masses = Hits$id.masses[ind.M]
    for (k in 1:nrow(pr)) {
        pr[k, ] <- pr[k, ]/sum(pr[k, ])
    }
    out <- list(Priors = pr, id.masses = id.masses, all.formulas = all.formulas)

    ## filter masses and formulas

    return(out)

}
