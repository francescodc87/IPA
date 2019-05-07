#' @title Function used to update your database
#'
#' @description
#' This function updates the database according to new information
#'
#' @param DB The database used in a dataframe format. It has to be organized in the same fashion of the one included with the package
#' @param all_adducts_matrix_POS A matrix contains information about all possible adducts.
#' Columns are: KEGG.id (e.g. 'C00002'), adduct (e.g. 'M+H'), RT (previously known
#' retention time for the compound, NA if unknown), formula (e.g. C10H17N5O13P3),
#' theoretical mz (e.g. 508.003570124), charge (e.g. 1) mono ('mono'
#' if monoisotopical form, 'iso' otherwise)
#' @param all_adducts_matrix_NEG A matrix contains information about all possible adducts.
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

"update_database" <- function(DB, all_adducts_matrix_POS, all_adducts_matrix_NEG, update, adducts, isotopes) {
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
                  all_adducts_matrix_POS <- rbind(all_adducts_matrix_POS, c(as.vector(tmp["KEGG.id"]), add, as.vector(tmp["RT"]), out, "mono"))
                }
            }
            ## NEG
            tmp.adds <- unlist(strsplit(tmp["NEG.adducts"], split = ";"))
            for (add in tmp.adds) {
                out <- Adducts.mz.formula(tmp.formula, add, adducts, isotopes)
                if (!is.null(out)) {
                  all_adducts_matrix_NEG <- rbind(all_adducts_matrix_NEG, c(as.vector(tmp["KEGG.id"]), add, as.vector(tmp["RT"]), out, "mono"))
                }
            }

        } else {
            for (field in colnames(update)) {
                DB[ind, field] <- update[k, field]
            }

            ### compute adducts
            tmp.formula <- as.vector(update[k, "Formula"])
            ## POS
            tmp.adds <- unlist(strsplit(update[k, "POS.adducts"], split = ";"))
            out.pos <-NULL
            for (add in tmp.adds) {
                out <- Adducts.mz.formula(tmp.formula, add, adducts, isotopes)
                if (!is.null(out)) {
                  out.pos <- rbind(out.pos, c(as.vector(update[k, "KEGG.id"]), add, as.vector(update[k, "RT"]), out, "mono"))
                }
            }
            ## NEG
            tmp.adds <- unlist(strsplit(update[k, "NEG.adducts"], split = ";"))
            out.neg <- NULL
            for (add in tmp.adds) {
                out <- Adducts.mz.formula(form = tmp.formula, add, adducts, isotopes)
                if (!is.null(out)) {
                  out.neg <- rbind(out.neg, c(as.vector(update[k, "KEGG.id"]), as.vector(update[k, "RT"]), add, out, "mono"))
                }
            }

            ### remove old adducts
            ind.tmp <- which(all_adducts_matrix_POS[, 1] == update[k, 1])
            if(ind.tmp[length(ind.tmp)]<nrow(all_adducts_matrix_POS)){
              all_adducts_matrix_POS <- rbind(all_adducts_matrix_POS[1:(ind.tmp[1]-1), ],out.pos,all_adducts_matrix_POS[(ind.tmp[length(ind.tmp)]+1):nrow(all_adducts_matrix_POS), ])
            }else{
              all_adducts_matrix_POS <- rbind(all_adducts_matrix_POS[1:(ind.tmp[1]-1), ],out.pos)
            }
            ind.tmp <- which(all_adducts_matrix_NEG[, 1] == update[k, 1])
            if(ind.tmp[length(ind.tmp)]<nrow(all_adducts_matrix_NEG)){
              all_adducts_matrix_NEG <- rbind(all_adducts_matrix_NEG[1:(ind.tmp[1]-1), ],out.neg,all_adducts_matrix_NEG[(ind.tmp[length(ind.tmp)]+1):nrow(all_adducts_matrix_NEG), ])
            }else{
              all_adducts_matrix_NEG <- rbind(all_adducts_matrix_NEG[1:(ind.tmp[1]-1), ],out.neg)
            }

        }
    }

    out <- list(DB = DB, all_adducts_matrix_POS = all_adducts_matrix_POS, all_adducts_matrix_NEG = all_adducts_matrix_NEG)
    return(out)

}

