"formula_isopattern" <- function(isotable, isotopes) {
    NM <- dim(isotable)
    rnames <- rep(NA, NM[1])
    new.names <- NULL
    for (i in 1:NM[1]) {
        tmp <- isotable[i, 3:NM[2]]
        tmp.elements <- names(tmp)
        tmp.elements.char <- (str_extract(tmp.elements, "[aA-zZ]+"))
        tmp.elements.num <- (str_extract(tmp.elements, "[0-9]+"))
        tmp.element.char.uni <- unique(tmp.elements.char)
        # hill system order
        if (any(tmp.element.char.uni == "C")) {
            ind.C <- which(tmp.element.char.uni == "C")
            ind.H <- which(tmp.element.char.uni == "H")
            tmp.element.char.uni <- c(tmp.element.char.uni[c(ind.C, ind.H)], sort(tmp.element.char.uni[!1:length(tmp.element.char.uni) %in% c(ind.C, ind.H)]))
        } else {
            tmp.element.char.uni <- sort(tmp.element.char.uni)
        }
        new_order <- NULL
        for (j in 1:length(tmp.element.char.uni)) {
            ind.tmp <- which(tmp.elements.char == tmp.element.char.uni[j])
            tmp.elements[ind.tmp[1]] <- tmp.elements.char[ind.tmp[1]]
            ind.tmp2 <- ind.tmp
            if (length(ind.tmp) > 1) {
                ind.tmp <- ind.tmp[2:length(ind.tmp)]
                tmp.elements[ind.tmp] <- paste("[", tmp.elements.num[ind.tmp], "]", tmp.elements.char[ind.tmp], sep = "")
            }
            new_order <- c(new_order, ind.tmp2)
        }
        tmp <- tmp[new_order]
        tmp.elements <- tmp.elements[new_order]
        tmp.elements.char <- tmp.elements.char[new_order]
        tmp.elements.num <- tmp.elements.num[new_order]
        keep <- which(tmp != 0)
        tmp <- paste(tmp.elements[keep], tmp[keep], sep = "")
        tmp <- paste(tmp, collapse = "")
        new.names <- c(new.names, tmp)
    }
    rownames(isotable) <- new.names
    return(isotable)
}
