#' @title Building Isotopes connectivity matrix
#'
#' @description
#' This functions created the Isotopes connectivity matrix needed for the posterior estimation
#'
#' @param Prior The output of compute.Priors() function
#' @param DB The database used in a dataframe format. It has to be organized in the same fashion of the one included with the package
#' @param ratios A logical value indicating if the matrix should contain the expected intensity ratios (default TRUE)
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 500)
#'
#' @return A matrix encoding the isotopes connections
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso build.add.connectivity.matrix build.bio.connectivity.matrix
#'
#' @import Matrix
#' @import enviPat
#' @import stringr
#' @importFrom utils combn
#' @export

"build.iso.connectivity.matrix" <- function(Prior, DB, ratios = TRUE, v = TRUE, IT = 500) {
    cat("Finding isotopes connections... \n")
    Iso <- Matrix(0, nrow(Prior$all.formulas), nrow(Prior$all.formulas))
    monos <- which(Prior$all.formulas[, 7] == "mono")
    N <- length(monos)

    if (ratios) {
        for (k in 1:N) {
            ind <- which(Prior$all.formulas[, 1] == Prior$all.formulas[monos[k], 1] & Prior$all.formulas[, 2] == Prior$all.formulas[monos[k], 2])
            if (length(ind) > 1) {
                abundances <- as.numeric(Prior$all.formulas[ind, 8])
                for (j in 1:(length(ind) - 1)) {
                  Iso[ind[j], ind[j + 1]] <- abundances[j]/abundances[j + 1]
                  Iso[ind[j + 1], ind[j]] <- abundances[j + 1]/abundances[j]
                }

            }
            if (v) {
                if (k%%IT == 0) {
                  # Print on the screen some message
                  cat(paste0(round(k/N, 1) * 100, "%", "\n"))
                }
            }

        }
    } else {
        for (k in 1:N) {
            ind <- which(Prior$all.formulas[, 1] == Prior$all.formulas[monos[k], 1] & Prior$all.formulas[, 2] == Prior$all.formulas[monos[k], 2])
            if (length(ind) > 1) {
                for (j in 1:(length(ind) - 1)) {
                  Iso[ind[j], ind[j + 1]] <- 1
                  Iso[ind[j + 1], ind[j]] <- 1
                }

            }
            if (v) {
                if (k%%IT == 0) {
                  # Print on the screen some message
                  cat(paste0(round(k/N, 1) * 100, "%", "\n"))
                }
            }

        }
    }

    return(Iso)

}
