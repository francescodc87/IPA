#' @title Finding putative annotations
#'
#' @description
#' This functions takes find all the putative annotations based
#'
#' @param adducts.matrix A matrix contains information about all possible adducts.
#' Columns are: KEGG.id (e.g. 'C00002'), adduct (e.g. 'M+H'), RT (previously known
#' retention time for the compound, NA if unknown), formula (e.g. C10H17N5O13P3),
#' theoretical mz (e.g. 508.003570124), charge (e.g. -1) mono ('mono'
#' if monoisotopical form, 'iso' otherwise)
#' @param dataset A matrix containing the measured data, organized in 3 colums: mz, RT and Int
#' @param ppm.thr A numerical value indicating the maximum accuracy value to be considered
#' @param RTwin1 A numerical value indicating the maximum difference allowed
#' between measured RT and previously known retention time for the compound
#' @param RTwin2 A numerical value indicating the maximum difference allowed
#' between measured RT of different isotopes (should belower than RTwin1)
#' @param isotopes A matrix containing infomation about isotopes
#' @param iso.threshold A numerical value indicating the probability below which
#' isotope peaks can be omitted
#' @param corr.matrix A matrix containing the correlation values between peaks
#' @param corr.thr A numerical value expressing the treshold used to consider the correlations significant
#' @param relation.id A vector containg class labels of the previously grouped peaks
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 120)
#'
#' @return A list containing the matrix of the posterior probabilities, the id.masses vector, the all.formulas dataframe and
#' the allsampcomp matrix containing all the assignments for each iteration of the Gibbs sampler
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

"find.hits" <- function(adducts.matrix, dataset, ppm.thr, RTwin1, RTwin2, isotopes, iso.threshold = 1, corr.matrix = NULL, corr.thr = 0.75, relation.id = NULL, v = T, IT = 500) {
    cat("Finding hits... \n")
    all.adducts.matrix.hits <- matrix(NA, 0, 9)
    colnames(all.adducts.matrix.hits) <- c("KEGG.id", "adduct", "RT", "Formula", "theoretical mz", "charge", "mono", "abundance", "mz.id")
    for (k in 1:nrow(dataset)) {
        mz <- as.numeric(dataset[k, 1])
        RT <- as.numeric(dataset[k, 2])
        deltaM <- ppm.thr * mz * 1e-06
        ind <- which(abs(as.numeric(adducts.matrix[, 5]) - mz) <= deltaM)
        if (length(ind) > 0) {
            for (h in 1:length(ind)) {
                ### for each hit I will compute the isotopes and look for hits!
                tmp <- t(c(adducts.matrix[ind[h], ], abundance = 100, mz.id = k))
                if (is.na(adducts.matrix[ind[h], 3])) {
                  RTdiff1 <- 0
                } else {
                  RTdiff1 <- abs(adducts.matrix[ind[h], 3] - RT)
                }
                if (RTdiff1 <= RTwin1) {
                  isotable <- isopattern(isotopes = isotopes, adducts.matrix[ind[h], 4], charge = as.numeric(adducts.matrix[ind[h], 6]), verbose = FALSE, threshold = iso.threshold)[[1]]
                  isotable <- formula_isopattern(isotable, isotopes)
                  if (nrow(isotable) > 1) {
                    isotable <- isotable[order(isotable[, 2], decreasing = T), ]
                  }
                  ### find hits
                  flag.iso = TRUE
                  count = 1
                  while (flag.iso & nrow(isotable) > 1 & count < nrow(isotable)) {
                    count = count + 1
                    deltaM.iso <- ppm.thr * (isotable[count, 1]) * 1e-06

                    if (is.null(corr.matrix) & is.null(relation.id)) {
                      ind.iso <- which(abs(as.numeric(dataset[, 1]) - isotable[count, 1]) <= deltaM.iso & abs(as.numeric(dataset[, 2]) - RT) <= RTwin2)
                    } else if (is.null(relation.id)) {
                      ind.iso <- which(abs(as.numeric(dataset[, 1]) - isotable[count, 1]) <= deltaM.iso & as.numeric(corr.matrix[k, ]) >= corr.thr)
                    } else {
                      ind.iso <- which(abs(as.numeric(dataset[, 1]) - isotable[count, 1]) <= deltaM.iso & relation.id == relation.id[k])
                    }

                    if (length(ind.iso) == 0) {
                      flag.iso <- FALSE
                    } else {
                      for (i in 1:length(ind.iso)) {
                        tmp.iso <- c(adducts.matrix[ind[h], 1:3], rownames(isotable)[count], isotable[count, 1], adducts.matrix[ind[h], 6], "iso", isotable[count, 2], mz.id = ind.iso[i])
                        tmp <- rbind(tmp, tmp.iso)
                      }


                    }
                  }
                }
                all.adducts.matrix.hits <- rbind(all.adducts.matrix.hits, tmp)
            }


        }
        if (v) {
            if (k%%IT == 0) {
                # Print on the screen some message
                cat(paste0(round(k/nrow(dataset), 1) * 100, "%", "\n"))
            }
        }

    }

    rownames(all.adducts.matrix.hits) <- NULL
    cat("Mapping results... \n")
    id.masses <- unique(all.adducts.matrix.hits[, 9])
    all.formulas <- unique(all.adducts.matrix.hits[, 1:8], MARGIN = 1)
    hit.matrix <- Matrix(0, nrow = length(id.masses), ncol = nrow(all.formulas))
    for (k in 1:length(id.masses)) {
        ind <- which(all.adducts.matrix.hits[, 9] == id.masses[k])
        tmp <- unique(all.adducts.matrix.hits[ind, 1:8])
        ind <- which(duplicated(rbind(all.formulas, tmp), fromLast = T))
        hit.matrix[k, ind] <- 1
        if (v) {
            if (k%%IT == 0) {
                # Print on the screen some message
                cat(paste0(round(k/length(id.masses), 1) * 100, "%", "\n"))
            }
        }

    }

    out <- list(id.masses = id.masses, all.formulas = all.formulas, hit.matrix = hit.matrix)
    return(out)
}
