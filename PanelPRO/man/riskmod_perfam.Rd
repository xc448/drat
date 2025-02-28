% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/ModPenetrance.R
\name{riskmod_perfam}
\alias{riskmod_perfam}
\title{Modify cancer penetrances based on hazard ratios or relative risks for interventions}
\usage{
riskmod_perfam(raw_famp, fam_distinct, riskmod)
}
\arguments{
\item{raw_famp}{A numeric array of cancer penetrances for individuals in
\code{fam_distinct}, created inside of the \code{calcCancerPenetrance} function.}

\item{fam_distinct}{A checked pedigree data frame returned by
\code{\link{checkFam}}.}

\item{riskmod}{Risk modifier information contained in the \code{riskmod}
component in the database returned by \code{\link{buildDatabase}}.}
}
\value{
A numeric array of modified cancer penetrances with the same
dimensions as \code{raw_famp}.
}
\description{
Modify cancer penetrances based on hazard ratios or relative risks for interventions
}
\details{
Risk modifiers are assumed to be independent.
}
\seealso{
Other modifypen: 
\code{\link{postIntervHR}()},
\code{\link{postIntervRR}()}
}
\concept{modifypen}
