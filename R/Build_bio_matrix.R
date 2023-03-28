#' @title Building Biotransformations connectivity matrix
#'
#' @description
#' This functions created the Biotransformations connectivity matrix needed for the posterior estimation
#'
#' @param Prior The output of compute.Priors() function
#' @param DB The database used in a dataframe format. It has to be organized in the same fashion of the one included with the package
#' @param ionisation A character indicating the ionisation mode of the experiment, either 'positive' or 'negative'
#' @param connection.type if equal to "reactions" (default) considers the reactions reported in the database. If it is equal to "formulas", the function consider the list of possible reactions reported in connections
#' @param connections if connection.type=="fomulas", this vector contains a list of formulas that are likely to be added or substracted by an enzymatic reaction. If it is equal to NULL the default list of formulas is used
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

"build.bio.connectivity.matrix" <- function(Prior, DB, ionisation, connection.type = "reactions", connections = NULL,v = TRUE, IT = 500) {
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

    reduced.DB <- DB[which(DB[, 1] %in% Prior$all.formulas[ind.mains, 1]), ]



    if(connection.type == "reactions"){
      cat("Finding biochemical connections based in database reactions... \n")
      all.reactions <- unique(unlist(strsplit(reduced.DB[, 6], split = "\\;| ")))
      all.reactions <- all.reactions[!is.na(all.reactions)]
      reactions.list <- strsplit(reduced.DB[, 6], split = "\\;| ")
      React.map <- matrix(0, length(ind.mains), length(all.reactions))
      for (k in 1:length(ind.mains)) {
        ind <- which(reduced.DB[, 1] == Prior$all.formulas[ind.mains[k], 1])
        if(length(ind)>0){
            tmp.reacts <- which(all.reactions %in% reactions.list[[ind]])
            React.map[k, tmp.reacts] <- 1
        }  
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

    }else if(connection.type == "formulas"){

      ### reorder reduced DB
      ordine <- NULL
      for(xx in 1:length(ind.mains)){
        ordine <- c(ordine,which(reduced.DB[,1]==Prior$all.formulas[ind.mains[xx],1]))
      }
      reduced.DB<- reduced.DB[ordine,]

      cat("Finding biochemical connections based on formulas... \n")
      if(is.null(connections)){
        connections <- c("C3H5NO", "C6H12N4O", "C4H6N2O2", "C4H5NO3", "C3H5NOS", "C6H10N2O3S2", "C5H7NO3",
                         "C5H8N2O2","C2H3NO","C6H7N3O","C6H11NO","C6H11NO","C6H12N2O","C5H9NOS","C9H9NO",
                         "C5H7NO","C3H5NO2","C4H7NO2","C11H10N2O","C9H9NO2","C5H9NO","C4H4O2","C3H5O",
                         "C10H12N5O6P","C10H15N2O3S","C10H14N2O2S","CH2ON","C21H34N7O16P3S","C21H33N7O15P3S",
                         "C10H15N3O5S","C5H7","C3H2O3","C16H30O","C8H8NO5P","CH3N2O","C5H4N5","C10H11N5O3",
                         "C10H13N5O9P2","C10H12N5O6P","C9H13N3O10P2","C9H12N3O7P","C4H4N3O","C10H13N5O10P2",
                         "C10H12N5O7P","C5H4N5O","C10H11N5O4","C10H14N2O10P2","C10H12N2O4","C5H5N2O2","C10H13N2O7P",
                         "C9H12N2O11P2","C9H11N2O8P","C4H3N2O2","C9H10N2O5","C2H3O2","C2H2O","C2H2","CO2","CHO2","H2O",
                         "H3O6P2","C2H4","CO","C2O2","H2","O","P","C2H2O","CH2","HPO3","NH2","PP","NH","SO3","N",
                         "C6H10O5","C6H10O6","C5H8O4","C12H20O11","C6H11O8P","C6H8O6","C6H10O5","C18H30O15") ### put the list I know alreadr
      }
      data("isotopes")
      connections <- check_chemform(isotopes = isotopes, chemforms = connections)
      connections.mw <- connections[,3]
      for(cc in 1:length(ind.mains)){
        tmp <- cbind(subform(reduced.DB[,4], reduced.DB[cc,4]),check_ded(reduced.DB[,4], reduced.DB[cc,4]))
        mw.diff <- rep(NA,nrow(tmp))
        mw.diff[which(tmp[,2]==FALSE)] <- check_chemform(isotopes = isotopes, chemforms = tmp[which(tmp[,2]==FALSE),1])[,3]
        tmp <- cbind(tmp, mw.diff)
        idx.sub <- which(tmp[,2]==FALSE & tmp[,1]!="NANA" & tmp[,3]%in%connections.mw)
        ccc <- ind.mains[cc]
        idx.sub <- ind.mains[idx.sub]
        if(length(idx.sub)>0){
          Bio[cbind(ccc,idx.sub)] <- 1
          Bio[cbind(idx.sub,ccc)] <- 1
        }
        if (v) {
          if (cc%%IT == 0) {
            # Print on the screen some message
            cat(paste0(round((cc * 100)/nrow(reduced.DB), 1), "%", "\n"))
          }
        }

      }



    }else{
      cat("\n wrong connection type")
      stop()
    }

    return(Bio)

}
