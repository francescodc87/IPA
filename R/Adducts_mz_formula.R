"Adducts.mz.formula" <- function(form, add, adducts, isotopes) {
    flag <- TRUE
    ### first check if the formula is in the correct format
    form <- check_chemform(isotopes, form)[1, 2]

    # First, multiply the chemical formula of the molecule by the times it appears in the final adduct; multiform
    ind <- which(adducts[, 1] == add)
    form <- multiform(form, adducts[ind, 4])

    # Second, add the chemical formula of any adduct to that of the molecule; mergeform if there is something to ADD
    if (adducts[ind, 7] != "FALSE") {
        form <- mergeform(form, adducts[ind, 7])
    }
    # Third, subtract the chemical formula of any deduct form that of the molecule; check_ded & subform. if there is something to subtract
    if (adducts[ind, 8] != "FALSE") {
        if (check_ded(form, adducts[ind, 8]) == "FALSE") {
            form <- subform(form, adducts[ind, 8])
        } else {
            flag <- FALSE
        }
    }
    if (flag) {
        mz <- check_chemform(isotopes, form)[1, 3]/abs(adducts[ind, 3])
        out <- c(form, mz, as.numeric(adducts[ind, 3]))
        names(out) <- c("formula", "m/z", "charge")
    } else {
        out <- NULL
    }

    return(out)
}
