#' @title Computing Posterior probabilities associated with putative annotations
#'
#' @description
#' This functions takes as input the output of the compute.Priors() function and estimates the
#' posterior probabilities for the putative annotations
#'
#' @param P  The output of compute.Priors() function
#' @param Add A binary matrix containing all the possible adducts connections
#' @param Iso A matrix containg the expected intensity ratios between the considered isotopes
#' @param Bio A binary matrix containing the the possible biochemical connections
#' @param Int A numerical vector containing the measure intensities
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
#' @param Corr.matrix A matrix containing the correlation values between peaks
#' @param corr.thr A numerical value expressing the treshold used to consider the correlations significant
#' @param RT A numerical vector containin the measured retention times
#' @param RT.win A numerical value expressing the maximum Retention time difference allowed between adducts/isotopes
#' @param rel.id A vector containg class labels of the previously grouped peaks
#' @param v A logical value indicating if the progress will be shown (default TRUE)
#' @param IT A number inticating after how many iteration an update should be shown (default 120)
#'
#' @return A list containing the matrix of the posterior probabilities, the id.masses vector, the all.formulas dataframe and
#' the allsampcomp matrix containing all the assignments for each iteration of the Gibbs sampler
#'
#' @author Francesco Del Carratore \email{francescodc87@@gmail.com}
#'
#' @seealso find.Hits compute.Priors
#'
#' @import Matrix
#' @import enviPat
#' @import stringr
#' @importFrom utils combn
#' @export


"IPAposteriors" <- function(P, Add = NULL, Iso = NULL, Bio = NULL, Int = NULL, delta.bio = 1, delta.add = 1, delta.iso = 1, allsampcomp = NULL, ratio.toll = 0.8, allsamp = FALSE,
    no.its = 1100, burn = 100, Corr.matrix = NULL, corr.thr = 0.8, RT = NULL, RT.win = 3, rel.id = NULL, v = TRUE, IT = 500) {

    if (is.null(Add) & is.null(Iso) & is.null(Bio)) {
        cat("\n Insert at least one connectivity matrix (Add, Iso and Bio)")
        stop()
    }

    if (!is.null(Iso) & is.null(Int)) {
        cat("\n Insert intensities to use Isotopes connections")
        stop()
    }

    ### ONLY BIO - - - IGNORING ALL POSSIBLE FILTERS
    if (is.null(Add) & is.null(Iso) & !is.null(Bio)) {
        cat("\n Computing Posterior probabilities based on Biotransformations only")
        out <- IPA.sampler.One_Binary_Matrix(P = P$Priors, Bio = Bio, no.its = no.its, burn = burn, delta = delta.bio, allsamp = allsamp, allsampcomp = allsampcomp, v = v,
            IT = IT)
    }

    ### ONLY ISO - - - selecting the best filter
    if (is.null(Add) & !is.null(Iso) & is.null(Bio)) {
        cat("\n Computing Posterior probabilities based on Isotopes only")
        if (!is.null(rel.id)) {
            cat("\n Considering relation ids \n")
            out <- IPA.sampler.One_Non_Binary_Matrix_relId(P = P$Priors, Iso = Iso, Int = Int, rel.id = rel.id, ratio.toll = ratio.toll, no.its = no.its, burn = burn, delta = delta.iso,
                allsamp = allsamp, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(Corr.matrix)) {
            cat("\n Considering correlation matrix \n")
            out <- IPA.sampler.One_Non_Binary_Matrix_CorrMat(P = P$Priors, Iso = Iso, Int = Int, Corr.matrix = Corr.matrix, corr.thr = corr.thr, ratio.toll = ratio.toll,
                no.its = no.its, burn = burn, delta = delta.iso, allsamp = allsamp, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(RT)) {
            cat("\n Considering Retention Times \n")
            out <- IPA.sampler.One_Non_Binary_Matrix_RT(P = P$Priors, Iso = Iso, Int = Int, RT = RT, RT.win = RT.win, ratio.toll = ratio.toll, no.its = no.its, burn = burn,
                delta = delta.iso, allsamp = allsamp, allsampcomp = allsampcomp, v = v, IT = IT)

        } else {
            cat("\n ERROR: please input relation ids or correlations matrix or retention times")
            stop()
        }

    }

    ### ONLY ADD - - - selecting the best filter
    if (!is.null(Add) & is.null(Iso) & is.null(Bio)) {
        cat("\n Computing Posterior probabilities based on Adducts only")
        if (!is.null(rel.id)) {
            cat("\n Considering relation ids \n")
            out <- IPA.sampler.One_Binary_Matrix_relId(P = P$Priors, Add = Add, rel.id = rel.id, no.its = no.its, burn = burn, delta = delta.add, allsamp = allsamp, allsampcomp = allsampcomp,
                v = v, IT = IT)
        } else if (!is.null(Corr.matrix)) {
            cat("\n Considering correlation matrix \n")
            out <- IPA.sampler.One_Binary_Matrix_CorrMat(P = P$Priors, Add = Add, Corr.matrix = Corr.matrix, corr.thr = corr.thr, no.its = no.its, burn = burn, delta = delta.add,
                allsamp = allsamp, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(RT)) {
            cat("\n Considering Retention Times \n")
            out <- IPA.sampler.One_Binary_Matrix_RT(P = P$Priors, Add = Add, RT = RT, RT.win = RT.win, no.its = no.its, burn = burn, delta = delta.add, allsamp = allsamp,
                allsampcomp = allsampcomp, v = v, IT = IT)

        } else {
            cat("\n ERROR: please input relation ids or correlations matrix or retention times")
            stop()
        }
    }


    ### ONLY ISO and ADD - - - selecting the best filter
    if (!is.null(Add) & !is.null(Iso) & is.null(Bio)) {
        cat("\n Computing Posterior probabilities based on Adducts and isotpes")
        if (!is.null(rel.id)) {
            cat("\n Considering relation ids \n")
            out <- IPA.sampler.Iso.Add.relID(P = P$Priors, Iso = Iso, Add = Add, Int = Int, rel.id = rel.id, no.its = no.its, burn = burn, delta.iso = delta.iso, delta.add = delta.add,
                allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(Corr.matrix)) {
            cat("\n Considering correlation matrix \n")
            out <- IPA.sampler.Iso.Add.CorrMat(P = P$Priors, Iso = Iso, Add = Add, Int = Int, Corr.matrix = Corr.matrix, corr.thr = corr.thr, no.its = no.its, burn = burn,
                delta.iso = delta.iso, delta.add = delta.add, allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(RT)) {
            cat("\n Considering Retention Times \n")
            out <- IPA.sampler.Iso.Add.RT(P = P$Priors, Iso = Iso, Add = Add, Int = Int, RT = RT, RT.win = RT.win, no.its = no.its, burn = burn, delta.iso = delta.iso, delta.add = delta.add,
                allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp, v = v, IT = IT)

        } else {
            cat("\n ERROR: please input relation ids or correlations matrix or retention times")
            stop()
        }
    }

    ### ONLY ISO and BIO - - - selecting the best filter
    if (is.null(Add) & !is.null(Iso) & !is.null(Bio)) {
        cat("\n Computing Posterior probabilities based on Isotopes and Biotransformations")
        if (!is.null(rel.id)) {
            cat("\n Considering relation ids \n")
            out <- IPA.sampler.Iso.Bio.relID(P = P$Priors, Iso = Iso, Bio = Bio, Int = Int, rel.id = rel.id, no.its = no.its, burn = burn, delta.iso = delta.iso, delta.bio = delta.bio,
                allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(Corr.matrix)) {
            cat("\n Considering correlation matrix \n")
            out <- IPA.sampler.Iso.Bio.CorrMat(P = P$Priors, Iso = Iso, Bio = Bio, Int = Int, Corr.matrix = Corr.matrix, corr.thr = corr.thr, no.its = no.its, burn = burn,
                delta.iso = delta.iso, delta.bio = delta.bio, allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(RT)) {
            cat("\n Considering Retention Times \n")
            out <- IPA.sampler.Iso.Bio.RT(P = P$Priors, Iso = Iso, Bio = Bio, Int = Int, RT = RT, RT.win = RT.win, no.its = no.its, burn = burn, delta.iso = delta.iso, delta.bio = delta.bio,
                allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp, v = v, IT = IT)

        } else {
            cat("\n ERROR: please input relation ids or correlations matrix or retention times")
            stop()
        }
    }


    ### ONLY ADD and BIO - - - selecting the best filter
    if (!is.null(Add) & is.null(Iso) & !is.null(Bio)) {
        cat("\n Computing Posterior probabilities based on Adducts and Biotransformations")
        if (!is.null(rel.id)) {
            cat("\n Considering relation ids \n")
            out <- IPA.sampler.Two_Binary_Matrix_relId_Add_Bio(P = P$Priors, Add = Add, Bio = Bio, rel.id = rel.id, no.its = no.its, burn = burn, delta.add = delta.add,
                delta.bio = delta.bio, allsamp = allsamp, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(Corr.matrix)) {
            cat("\n Considering correlation matrix \n")
            out <- IPA.sampler.Two_Binary_Matrix_CorrMat_Add_Bio(P = P$Priors, Add = Add, Bio = Bio, Corr.matrix = Corr.matrix, corr.thr = corr.thr, no.its = no.its, burn = burn,
                delta.add = delta.add, delta.bio = delta.bio, allsamp = allsamp, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(RT)) {
            cat("\n Considering Retention Times \n")
            out <- IPA.sampler.Two_Binary_Matrix_RT_Add_Bio(P = P$Priors, Add = Add, Bio = Bio, RT = RT, RT.win = RT.win, no.its = no.its, burn = burn, delta.add = delta.add,
                delta.bio = delta.bio, allsamp = allsamp, allsampcomp = allsampcomp, v = v, IT = IT)

        } else {
            cat("\n ERROR: please input relation ids or correlations matrix or retention times")
            stop()
        }
    }


    ### ISO ADD and BIO - - - selecting the best filter
    if (!is.null(Add) & !is.null(Iso) & !is.null(Bio)) {
        cat("\n Computing Posterior probabilities based on Isotopes, Adducts and Biotransformations")
        if (!is.null(rel.id)) {
            cat("\n Considering relation ids \n")
            out <- IPA.sampler.Add.Iso.Bio.relID(P = P$Priors, Add = Add, Iso = Iso, Bio = Bio, Int = Int, rel.id = rel.id, no.its = no.its, burn = burn, delta.add = delta.add,
                delta.iso = delta.iso, delta.bio = delta.bio, allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp, v = v, IT = IT)

        } else if (!is.null(Corr.matrix)) {
            cat("\n Considering correlation matrix \n")
            out <- IPA.sampler.Add.Iso.Bio.CorrMat(P = P$Priors, Add = Add, Iso = Iso, Bio = Bio, Int = Int, Corr.matrix = Corr.matrix, corr.thr = corr.thr, no.its = no.its,
                burn = burn, delta.add = delta.add, delta.iso = delta.iso, delta.bio = delta.bio, allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp,
                v = v, IT = IT)

        } else if (!is.null(RT)) {
            cat("\n Considering Retention Times \n")
            out <- IPA.sampler.Add.Iso.Bio.RT(P = P$Priors, Add = Add, Iso = Iso, Bio = Bio, Int = Int, RT = RT, RT.win = RT.win, no.its = no.its, burn = burn, delta.add = delta.add,
                delta.iso = delta.iso, delta.bio = delta.bio, allsamp = allsamp, ratio.toll = ratio.toll, allsampcomp = allsampcomp, v = v, IT = IT)

        } else {
            cat("\n ERROR: please input relation ids or correlations matrix or retention times")
            stop()
        }
    }



    if (allsamp) {
        Post <- list(Post = out$Post, id.masses = P$id.masses, all.formulas = P$all.formulas, allsampcomp = out$allsampcomp)
    } else {
        Post <- list(Post = out, id.masses = P$id.masses, all.formulas = P$all.formulas)
    }
    return(Post)
}
