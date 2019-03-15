#' @title Function used to update your database
#'
#' @description
#' This function updates the database according to new information
#'
#' @param DB The database used in a dataframe format. It has to be organized in the same fashion of the one included with the package
#' @param all.adducts.matrix.POS A matrix contains information about all possible adducts.
#' Columns are: KEGG.id (e.g. 'C00002'), adduct (e.g. 'M+H'), RT (previously known
#' retention time for the compound, NA if unknown), formula (e.g. C10H17N5O13P3),
#' theoretical mz (e.g. 508.003570124), charge (e.g. 1) mono ('mono'
#' if monoisotopical form, 'iso' otherwise)
#' @param all.adducts.matrix.NEG A matrix contains information about all possible adducts.
#' Columns are: KEGG.id (e.g. 'C00002'), adduct (e.g. 'M-H'), RT (previously known
#' retention time for the compound, NA if unknown), formula (e.g. C10H17N5O13P3),
#' theoretical mz (e.g. 508.003570124), charge (e.g. -1) mono ('mono'
#' if monoisotopical form, 'iso' otherwise)
#' @param update A dataframe containing the new information
#' @param adducts A dataframe containing adducts information as the one included in the enviPat package
#' @param isotopes A dataframe containing isotopes information as the one included in the enviPat package
#'
#' @return Updated database
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#'
#' @import Matrix
#' @import enviPat
#' @import stringr
#' @importFrom utils combn
#' @export

"update_database" <- function(DB, all.adducts.matrix.POS, all.adducts.matrix.NEG, update, adducts, isotopes) {
    for (k in 1:nrow(update)) {
        ind <- which(DB[, 1] == update[k, 1])
        if (length(ind) == 0) {

            tmp <- rep(NA, ncol(DB))
            names(tmp) <- colnames(DB)
            for (field in colnames(update)) {
                tmp[field] <- update[k, field]
            }
            DB <- rbind(DB, tmp)

            ### compute adducts
            tmp.formula <- as.vector(tmp["Formula"])
            ## POS
            tmp.adds <- unlist(strsplit(tmp["POS.adducts"], split = ";"))
            for (add in tmp.adds) {
                out <- Adducts.mz.formula(tmp.formula, add, adducts, isotopes)
                if (!is.null(out)) {
                  all.adducts.matrix.POS <- rbind(all.adducts.matrix.POS, c(as.vector(tmp["KEGG.id"]), add, as.vector(tmp["RT"]), out, "mono"))
                }
            }
            ## NEG
            tmp.adds <- unlist(strsplit(tmp["NEG.adducts"], split = ";"))
            for (add in tmp.adds) {
                out <- Adducts.mz.formula(tmp.formula, add, adducts, isotopes)
                if (!is.null(out)) {
                  all.adducts.matrix.NEG <- rbind(all.adducts.matrix.NEG, c(as.vector(tmp["KEGG.id"]), add, as.vector(tmp["RT"]), out, "mono"))
                }
            }

        } else {
            for (field in colnames(update)) {
                DB[ind, field] <- update[k, field]
            }
            ### remove old adducts
            ind.tmp <- which(all.adducts.matrix.POS[, 1] != update[k, 1])
            all.adducts.matrix.POS <- all.adducts.matrix.POS[ind.tmp, ]
            ind.tmp <- which(all.adducts.matrix.NEG[, 1] != update[k, 1])
            all.adducts.matrix.NEG <- all.adducts.matrix.NEG[ind.tmp, ]
            ### compute adducts
            tmp.formula <- as.vector(update[k, "Formula"])
            ## POS
            tmp.adds <- unlist(strsplit(update[k, "POS.adducts"], split = ";"))
            for (add in tmp.adds) {
                out <- Adducts.mz.formula(tmp.formula, add, adducts, isotopes)
                if (!is.null(out)) {
                  all.adducts.matrix.POS <- rbind(all.adducts.matrix.POS, c(as.vector(update[k, "KEGG.id"]), add, as.vector(update[k, "RT"]), out, "mono"))
                }
            }
            ## NEG
            tmp.adds <- unlist(strsplit(update[k, "NEG.adducts"], split = ";"))
            for (add in tmp.adds) {
                out <- Adducts.mz.formula(tmp.formula, add, adducts, isotopes)
                if (!is.null(out)) {
                  all.adducts.matrix.NEG <- rbind(all.adducts.matrix.NEG, c(as.vector(update[k, "KEGG.id"]), as.vector(update[k, "RT"]), add, out, "mono"))
                }
            }

        }
    }

    out <- list(DB = DB, all.adducts.matrix.POS = all.adducts.matrix.POS, all.adducts.matrix.NEG = all.adducts.matrix.NEG)
    return(out)

}

