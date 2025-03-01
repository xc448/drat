% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/RcppExports.R
\name{.peelingParing}
\alias{.peelingParing}
\title{Peeling and Paring}
\usage{
.peelingParing(
  probandIdxs,
  idMatrix,
  mates,
  children,
  prevalence,
  likelihood,
  transition,
  antProb,
  postProb
)
}
\arguments{
\item{probandIdxs}{The idx numbers of the probands in the ID column}

\item{idMatrix}{The matrix representing the idx numbers of the ID, MotherID and FatherID}

\item{mates}{The `matrix' specifying the mates}

\item{children}{The `cube' specifying the children}

\item{prevalence}{The `matrix' specifying the prevalences}

\item{likelihood}{The `matrix' specifying the likelihoods (penetrances)}

\item{transition}{The `cube' specifying the transition probabilities from mother and father to child}

\item{antProb}{A placeholder to memoise the recursive function output of calAntProb}

\item{postProb}{A placeholder to memoise the recursive function output of calPostProb}
}
\value{
Matrix of the posterior probabilities of ID having certain genotypes
}
\description{
Calculates the posterior probabilities of genotypes of probands given inputs
}
\details{
Calls calAntProb and calPostProb recursively
NB: Require that NAs in antProb and postProb are set to 0
}
