#' @title Building Biotransformations connectivity matrix
#'
#' @description
#' This functions created the Biotransformations connectivity matrix needed for the posterior estimation
#'
#' @param Prior The output of compute.Priors() function
#' @param DB The database used in a dataframe format. It has to be organized in the same fashion of the one included with the package
#' @param ionisation A character indicating the ionisation mode of the experiment, either 'positive' or 'negative'
#' @param connection.type NOT USED YET!
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 500)
#'
#' @return A binary matrix encoding the biotransfiormations connections
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso build.add.connectivity.matrix build.iso.connectivity.matrix
#' @import Matrix
#' @import enviPat
#' @import stringr
#' @importFrom utils combn
#' @export

"build.bio.connectivity.matrix" <- function(Prior, DB, ionisation, connection.type = "reactions", v = TRUE, IT = 500) {
    Bio <- Matrix(0, nrow(Prior$all.formulas), nrow(Prior$all.formulas))
    IDs <- unique(Prior$all.formulas[, 1])
    cat("Finding main adducts... \n")
    ind.mains <- NULL
    for (k in 1:length(IDs)) {
        if (ionisation == "positive") {
            main <- DB[DB[, 1] == IDs[k], 12]
        } else {
            main <- DB[DB[, 1] == IDs[k], 14]
        }
        ind <- which(Prior$all.formulas[, 1] == IDs[k] & Prior$all.formulas[, 2] == main & Prior$all.formulas[, 7] == "mono")
        if (length(ind) < 1) {
            ind <- which(Prior$all.formulas[, 1] == IDs[k] & Prior$all.formulas[, 7] == "mono")[1]
        }
        ind.mains <- c(ind.mains, ind)
        if (v) {
            if (k%%IT == 0) {
                # Print on the screen some message
                cat(paste0(round((k * 100)/length(IDs), 1), "%", "\n"))
            }
        }
    }

    ind.mains <- ind.mains[!is.na(ind.mains)]
    ### now I have to map the reactions based on the kegg reactions or on the list of mass differences!
    cat("Finding biochemical connections... \n")
    reduced.DB <- DB[which(DB[, 1] %in% Prior$all.formulas[ind.mains, 1]), ]
    all.reactions <- unique(unlist(strsplit(reduced.DB[, 6], split = "\\;| ")))
    all.reactions <- all.reactions[!is.na(all.reactions)]
    reactions.list <- strsplit(reduced.DB[, 6], split = "\\;| ")
    React.map <- matrix(0, length(ind.mains), length(all.reactions))
    for (k in 1:length(ind.mains)) {
        ind <- which(reduced.DB[, 1] == Prior$all.formulas[ind.mains[k], 1])
        tmp.reacts <- which(all.reactions %in% reactions.list[[ind]])
        React.map[k, tmp.reacts] <- 1
    }

    for (k in 1:ncol(React.map)) {
        ind <- ind.mains[which(React.map[, k] == 1)]
        if (length(ind) == 2) {
            Bio[ind[1], ind[2]] <- 1
            Bio[ind[2], ind[1]] <- 1
        } else if (length(ind) > 2) {
            ind <- t(combn(ind, 2))
            Bio[ind] <- 1
            Bio[ind[, 2:1]] <- 1
        }
        if (v) {
            if (k%%IT == 0) {
                # Print on the screen some message
                cat(paste0(round((k * 100)/length(all.reactions), 1), "%", "\n"))
            }
        }
    }

    return(Bio)


}
