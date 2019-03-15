#' @title Building Adducts connectivity matrix
#'
#' @description
#' This functions created the Adducts connectivity matrix needed for the posterior estimation
#'
#' @param Prior The output of compute.Priors() function
#' @param DB The database used in a dataframe format. It has to be organized in the same fashion of the one included with the package
#' @param ionisation A character indicating the ionisation mode of the experiment, either 'positive' or 'negative'
#' @param fully.connected A logical value if TRUE all adducts will be considered connected with each other, if FALSE (default) all
#' the adducts are only connected to the main adduct
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 500)
#'
#' @return A binary matrix encoding the adducts connections
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso build.bio.connectivity.matrix build.iso.connectivity.matrix
#' @import Matrix
#' @import enviPat
#' @import stringr
#' @importFrom utils combn
#' @export

"build.add.connenctivity.matrix" <- function(Prior, DB, ionisation, fully.connected = FALSE, v = TRUE, IT = 500) {
    cat("Finding adducts connections... \n")
    Add <- Matrix(0, nrow(Prior$all.formulas), nrow(Prior$all.formulas))
    IDs <- unique(Prior$all.formulas[, 1])
    N <- length(IDs)
    if (fully.connected) {
        for (k in 1:N) {
            ind <- which(Prior$all.formulas[, 1] == IDs[k] & Prior$all.formulas[, 7] == "mono")
            if (length(ind) > 1) {
                ind <- combn(ind, 2)
                ind <- rbind(t(ind), t(ind[2:1, ]))
                Add[ind] <- 1
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
            if (ionisation == "positive") {
                main <- DB[DB[, 1] == IDs[k], 12]
            } else {
                main <- DB[DB[, 1] == IDs[k], 14]
            }
            ind2 <- which(Prior$all.formulas[, 1] == IDs[k] & Prior$all.formulas[, 7] == "mono" & Prior$all.formulas[, 2] == main)
            ind <- which(Prior$all.formulas[, 1] == IDs[k] & Prior$all.formulas[, 7] == "mono" & Prior$all.formulas[, 2] != main)
            if (length(ind2) == 1 & length(ind) > 0) {
                ind <- rbind(cbind(ind, ind2), cbind(ind2, ind))
                Add[ind] <- 1
            }

            if (v) {
                if (k%%IT == 0) {
                  # Print on the screen some message
                  cat(paste0(round(k/N, 1) * 100, "%", "\n"))
                }
            }

        }
    }

    return(Add)

}
