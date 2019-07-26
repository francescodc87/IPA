#' @title Single function for the whole IPA analysis
#'
#' @description
#' Wrapper summmarizing the whole IPA analysis in one single function. USE IT CAREFULLY!
#'
#'
#' @param dataset A matrix containing the measured data, organized in 3 colums: mz, RT and Int
#' @param adducts.matrix A matrix contains information about all possible adducts.
#' Columns are: KEGG.id (e.g. 'C00002'), adduct (e.g. 'M+H'), RT (previously known
#' retention time for the compound, NA if unknown), formula (e.g. C10H17N5O13P3),
#' theoretical mz (e.g. 508.003570124), charge (e.g. -1) mono ('mono'
#' if monoisotopical form, 'iso' otherwise)
#' @param DB The database used in a dataframe format. It has to be organized in the same fashion of the one included with the package
#' @param ionisation A character indicating the ionisation mode of the experiment, either 'positive' or 'negative'
#' @param ppm A number indicating the instrument accuracy to be considered
#' @param ppm.thr A numerical value indicating the maximum accuracy value to be considered
#' @param RTwin A numerical value indicating the maximum difference allowed
#' between measured RT of different isotopes
#' @param isotopes A matrix containing infomation about isotopes
#' @param iso A logical parameter. If true isotope connections are used
#' @param add A logical parameter. If true adducts connections are used
#' @param bio A logical parameter. If true biochemical connections are used
#' @param ints A logical parameter. If true intensities are used
#' @param rt A logical parameter. If true RTs are used for filtering connections
#' @param iso.threshold A numerical value indicating the probability below which
#' isotope peaks can be omitted
#' @param corr.matrix A matrix containing the correlation values between peaks
#' @param corr.thr A numerical value expressing the treshold used to consider the correlations significant
#' @param relation.id A vector containg class labels of the previously grouped peaks
#' @param pk A vector of length nrow(adducts.matrix). Values between 0 and 1 expressing initial confidence of the presence of each formula
#' @param delta.bio A number expressing the confidence on the information encoded in Bio
#' (smaller the value higher the confidence)
#' @param delta.add A number expressing the confidence on the information encoded in Add
#' (smaller the value higher the confidence)
#' @param delta.iso A number expressing the confidence on the information encoded in Iso
#' (smaller the value higher the confidence)
#' @param allsampcomp A matrix containing all the assignments for each iteration of the Gibbs sampler in a previus run
#' @param ratio.toll A numerical value expressing the minimum accepeted ratio between thereoretical
#' and observed intensity ratios between isotopes (default 0.8)
#' @param allsamp A logical value indicating if the output should contain the allsampcomp object
#' @param no.its A numerical value indicating the number of iterations the to be performed by the Gibbs sampler (default 1100)
#' @param burn A numerical value indicating the number of initial iterations to be ignored when computing the posterior
#' probabilities (default 100)
#' @param RT.win.post A numerical value expressing the maximum Retention time difference allowed between adducts/isotopes
#' @param IDs name of each feature in the initial dataset
#' @param unknown.ppm The ppm number to be assigned to unknown (default NA - No unknown is considered)
#' @param RT.pen A numerical value indicating the multiplicative factor used when the measured retention time is outside of the range reported in the database
#' @param pr.lim A number inidicating the lowest probability value allowed (default 1e-05)
#' @param fully.connected A logical value if TRUE all adducts will be considered connected with each other, if FALSE (default) all
#' the adducts are only connected to the main adduct
#' @param ratios A logical value indicating if the matrix should contain the expected intensity ratios (default TRUE)
#' @param connection.type if equal to "reactions" (default) considers the reactions reported in the database. If it is equal to "formulas", the function consider the list of possible reactions reported in connections
#' @param connections if connection.type=="fomulas", this vector contains a list of formulas that are likely to be added or substracted by an enzymatic reaction. If it is equal to NULL the default list of formulas is used
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 120)
#'
#'
#' @return A list containing the Hits, Prior, Post and parsed results
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


"simpleIPA" <- function(dataset, adducts.matrix, DB, ionisation, ppm, ppm.thr, RTwin,isotopes, iso=TRUE, add=TRUE, bio=TRUE, ints=TRUE, rt=TRUE, iso.threshold=1,
                        corr.matrix = NULL, corr.thr = 0.75, relation.id = NULL, pk=NULL, delta.bio = 1, delta.add = 1, delta.iso = 1, allsampcomp = NULL,
                        ratio.toll = 0.8, allsamp = FALSE, no.its = 1100, burn = 100, RT.win.post = 3, IDs=NULL,
                        unknown.ppm = NA, RT.pen=0.5, pr.lim = 1e-05, fully.connected = FALSE, ratios = TRUE, connection.type = "reactions", connections = NULL,
                        v = T, IT = 500) {

  cat("\n this wrapper should be used carefully, PLEASE CONSIDER AVOIDING IT" )

  Hits <- find.hits(adducts.matrix= adducts.matrix, dataset=dataset,
                    ppm.thr= ppm.thr, RTwin=RTwin,
                    isotopes=isotopes, iso.threshold = iso.threshold,
                    corr.matrix = corr.matrix, corr.thr = corr.thr,
                    relation.id = relation.id,
                    v = v, IT = IT)

  if(is.null(pk)){
    pk <- rep(1, nrow(Hits$all.formulas))
  }else{
    pk <- pk[Hits$id.masses]
  }
  Prior <- compute.Priors(Hits=Hits, dataset=dataset,
                          pk=rep(1,nrow(Hits$all.formulas)),
                          ppm=ppm,unknown.ppm = unknown.ppm,
                          RT.pen=RT.pen, v = v, IT = IT, pr.lim = pr.lim)


  ### building ADD matrix
  if(add){
    ADD <- build.add.connenctivity.matrix(Prior=Prior,  DB=DB,
                                          ionisation=ionisation, fully.connected=fully.connected,
                                          v = v, IT = IT)
  }else{
    ADD <-NULL
  }


  ### building ISO matrix
  if(iso){
    ISO <- build.iso.connectivity.matrix(Prior=Prior, DB=DB, ratios=ratios, v = v, IT = IT)
  }else{
    ISO <- NULL
  }

  ### building BIO matrix
  if(bio){
    BIO1 <- build.bio.connectivity.matrix(Prior=Prior, DB=DB,
                                          ionisation=ionisation,
                                          connection.type=connection.type,
                                          connections = connections,
                                          v = v, IT = IT)
  }else{
    BIO1 <- NULL
  }

  if(ints){
    Int = as.numeric(dataset[,3])
  }else{
    Int = NULL
  }
  if(rt){
    RT = as.numeric(dataset[,2])
  }else{
    RT = NULL
  }



  Post <- IPAposteriors(P=Prior,Iso = ISO, Add = ADD, Bio = BIO1, Int = Int,
                        delta.bio = delta.bio, delta.add = delta.add, delta.iso = delta.iso, allsampcomp = allsampcomp,
                        ratio.toll = ratio.toll, allsamp = allsamp,
                        no.its = no.its, burn = burn, Corr.matrix =corr.matrix, corr.thr = corr.thr, RT = NULL, RT.win = RT.win.post,
                        rel.id = relation.id, v = v, IT = IT)

  Final.res <- ParseIPAresults(Post=Post, Prior=Prior,dataset =dataset, DB = DB, IDs=IDs, v = v, IT = IT)

  out <- list(Hits=Hits, Prior= Prior, Post=Post, Final.res=Final.res)
  return(out)
}
