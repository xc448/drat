% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/calcFutureRisk.R
\name{fr.cbc}
\alias{fr.cbc}
\title{Calculate Future Risk of Contralateral Breast Cancer (CBC)}
\usage{
fr.cbc(pen.cbc.net, pen.d.crude, pen.bc1.net, age.bc1, age.current, return.NA)
}
\arguments{
\item{pen.cbc.net}{Matrix of net CBC penetrances where row names are
genotypes and column names are ages from 1 to MAXAGE.}

\item{pen.d.crude}{Matrix of crude death from non-cancer cause penetrances
where rows names are genotypes and column names are ages from 1 to MAXAGE.}

\item{pen.bc1.net}{1 row matrix or a vector of net 1st Breast Cancer (BC)
penetrances from age 1 to MAXAGE. Penetrances should be based on a weighted
average of possible genotypes.}

\item{age.bc1}{Age at which the proband developed 1st BC; numeric value.}

\item{age.current}{The proband's current age; numeric value.}

\item{return.NA}{logical which indicates if a matrix of all NAs should be
returned instead of calculated values. Useful those without a first BC.}
}
\value{
A matrix of future risks where the row names are genotypes and the column
names are ages 1 to MAXAGE.
}
\description{
Calculates the future risk for developing cancer in the second breast at a
given future age for a subject previously affected with breast cancer. This
calculation considers both survival from death from non-cancer causes
and survival from CBC.
}
