% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/docData.R
\docType{data}
\name{log.rr.bcrat}
\alias{log.rr.bcrat}
\title{Coefficients in BCRAT relative risk model}
\format{
Race-specific coefficients from Gail et al. (1989) for White women,
Gail et al. (2007) for Black women, Matsuno et al. (2011) for Asian women,
and Banegas et al. (2017) for Hispanic women, saved as a data frame with 6
rows and 7 columns:
\itemize{
\item \code{beta1}: Coefficient for number of benign breast biopsies.
\item \code{beta2}: Coefficient for age at menarche.
\item \code{beta3}: Coefficient for age at first live birth.
\item \code{beta4}: Coefficient for number of female first-degree relatives with
breast cancer.
\item \code{beta5}: Coefficient for interaction between number of benign breast
biopsies and age>=50.
\item \code{beta6}: Coefficient for interaction between age at first live birth and
number of female first-degree relatives with breast cancer.
\item \code{Race}: \code{"White"}, \code{"Black"}, \code{"Hispanic"}, \code{"FHispanic"}, \code{"Asian"}, or
\code{"Other"}; \code{"Hispanic"} corresponds to U.S.-born Hispanic, while
\code{"FHispanic"} corresponds to foreign-born Hispanic.
}
}
\usage{
log.rr.bcrat
}
\description{
Coefficients in BCRAT relative risk model
}
\references{
Gail MH, Brinton LA, Byar DP, Corle DK, Green SB, Schairer C, Mulvihill JJ.
Projecting individualized probabilities of developing breast cancer for
white females who are being examined annually. JNCI: Journal of the National
Cancer Institute. 1989 Dec 20;81(24):1879-86.

Gail MH, Costantino JP, Pee D, Bondy M, Newman L, Selvan M, Anderson GL,
Malone KE, Marchbanks PA, McCaskill-Stevens W, Norman SA. Projecting
individualized absolute invasive breast cancer risk in African American
women. Journal of the National Cancer Institute. 2007 Dec 5;99(23):1782-92.

Matsuno RK, Costantino JP, Ziegler RG, Anderson GL, Li H, Pee D, Gail MH.
Projecting individualized absolute invasive breast cancer risk in Asian and
Pacific Islander American women. Journal of the National Cancer Institute.
2011 Jun 22;103(12):951-61.

Banegas MP, John EM, Slattery ML, Gomez SL, Yu M, LaCroix AZ, Pee D,
Chlebowski RT, Hines LM, Thompson CA, Gail MH. Projecting individualized
absolute invasive breast cancer risk in US Hispanic women. Journal of the
National Cancer Institute. 2017 Feb 1;109(2):djw215.
}
\keyword{datasets}
