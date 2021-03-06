% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/IPA_Posteriors.R
\name{IPAposteriors}
\alias{IPAposteriors}
\title{Computing Posterior probabilities associated with putative annotations}
\usage{
IPAposteriors(P, Add = NULL, Iso = NULL, Bio = NULL, Int = NULL,
  delta.bio = 1, delta.add = 1, delta.iso = 1, allsampcomp = NULL,
  ratio.toll = 0.8, allsamp = FALSE, no.its = 1100, burn = 100,
  Corr.matrix = NULL, corr.thr = 0.8, RT = NULL, RT.win = 3,
  rel.id = NULL, v = TRUE, IT = 500)
}
\arguments{
\item{P}{The output of compute.Priors() function}

\item{Add}{A binary matrix containing all the possible adducts connections}

\item{Iso}{A matrix containg the expected intensity ratios between the considered isotopes}

\item{Bio}{A binary matrix containing the the possible biochemical connections}

\item{Int}{A numerical vector containing the measure intensities}

\item{delta.bio}{A number expressing the confidence on the information encoded in Bio
(smaller the value higher the confidence)}

\item{delta.add}{A number expressing the confidence on the information encoded in Add
(smaller the value higher the confidence)}

\item{delta.iso}{A number expressing the confidence on the information encoded in Iso
(smaller the value higher the confidence)}

\item{allsampcomp}{A matrix containing all the assignments for each iteration of the Gibbs sampler in a previus run}

\item{ratio.toll}{A numerical value expressing the minimum accepeted ratio between thereoretical
and observed intensity ratios between isotopes (default 0.8)}

\item{allsamp}{A logical value indicating if the output should contain the allsampcomp object}

\item{no.its}{A numerical value indicating the number of iterations the to be performed by the Gibbs sampler (default 1100)}

\item{burn}{A numerical value indicating the number of initial iterations to be ignored when computing the posterior
probabilities (default 100)}

\item{Corr.matrix}{A matrix containing the correlation values between peaks}

\item{corr.thr}{A numerical value expressing the treshold used to consider the correlations significant}

\item{RT}{A numerical vector containin the measured retention times}

\item{RT.win}{A numerical value expressing the maximum Retention time difference allowed between adducts/isotopes}

\item{rel.id}{A vector containg class labels of the previously grouped peaks}

\item{v}{A logical value indicating if the progress will be shown (default TRUE)}

\item{IT}{A number inticating after how many iteration an update should be shown (default 120)}
}
\value{
A list containing the matrix of the posterior probabilities, the id.masses vector, the all.formulas dataframe and
the allsampcomp matrix containing all the assignments for each iteration of the Gibbs sampler
}
\description{
This functions takes as input the output of the compute.Priors() function and estimates the
posterior probabilities for the putative annotations
}
\seealso{
find.Hits compute.Priors
}
\author{
Francesco Del Carratore \email{francescodc87@gmail.com}
}
