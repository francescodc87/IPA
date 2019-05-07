#' @title Parse IPA resutls
#'
#' @description
#' This functions takes as input the output of the IPAposteriors() function and parses
#' the results in a more readable fashion
#'
#' @param Post  The output of IPAposteriors() function
#' @param Prior The output of compute.Priors() function
#' @param dataset The initial dataset
#' @param IDs name of each feature in the initial dataset
#' @param DB Database
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 120)
#'
#' @return A list containing...NOT SURE YET
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso find.Hits compute.Priors IPAposteriors
#'
#' @import Matrix
#' @import enviPat
#' @import stringr
#' @importFrom utils combn
#' @export


"ParseIPAresults" <- function(Post, Prior, dataset, DB, IDs=NULL, v = TRUE, IT = 500) {

  Final.res <- vector("list", length(Post$id.masses))
  if(!is.null(IDs)){
    IDs <- IDs[Post$id.masses]
  }
  #### I am here... prepare the table perfectly!!

  for(k in 1:length(Post$id.masses)){
    ind <- which(Prior$Priors[k,]>0)
    if(length(ind)==1){
      out <- t(c(Prior$all.formulas[ind,c(1,2,4,5,8)], Prior$Priors[k,ind], Post$Post[k,ind]))
    }else{
      out <- cbind(Prior$all.formulas[ind,c(1,2,4,5,8)], Prior$Priors[k,ind], Post$Post[k,ind])
    }
    nomi <- NULL
    ppms <- NULL
    for(j in 1:nrow(out)){
      idx <- which(DB[,1]==out[j,1])
      if(length(idx)==0){
        nomi<- c(nomi, NA)
      }else{
        nomi <- c(nomi, DB[idx,"Name"])
      }
      delta <- dataset[as.numeric(Prior$id.masses[k]),1]-as.numeric(out[j,4])
      ppms <- c(ppms, round((delta/as.numeric(out[j,4]))*1e6,4))
    }

    out <- cbind(out,nomi, ppms)
    out <- out[,c(1,8,2:4,9,5:7)]
    out <- matrix(out, ncol=9)
    colnames(out) <- c("KEGG id", "name",  "adduct", "Formula", "theoretical mz","ppm", "abundance","Prior", "Post")
    rownames(out) <- NULL
    Final.res[[k]] <- out
  }

  names(Final.res) <- IDs

  return(Final.res)

}
