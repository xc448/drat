% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/checkFamHelpers.R
\name{.checkMating}
\alias{.checkMating}
\title{Check a pedigree's mating structure for validity}
\usage{
.checkMating(fff, rel_l)
}
\arguments{
\item{fff}{A pedigree data frame.}

\item{rel_l}{The list of relationships in \code{ped} returned by a call to
\code{\link{get_relation_list}}.}
}
\value{
None
}
\description{
Raises errors for invalid mating patterns, i.e. intermarriages that create
"loops" in the family:
\itemize{
\item Mating between siblings or first cousins.
\item Mating between a child's paternal aunt/uncle and maternal aunt/uncle.
\item Mating between children and their parents, aunts, or uncles.
\item Mating between a person's relatives and the relatives of any of their
spouses.
}
}
\seealso{
Other check: 
\code{\link{.boundAffectionAges}()},
\code{\link{.checkGender}()},
\code{\link{.checkHeredity}()},
\code{\link{.checkTwins}()},
\code{\link{.findIDs}()},
\code{\link{.fix_gender_cancer}()},
\code{\link{.getOtherRels}()},
\code{\link{.getSpouses}()},
\code{\link{.pedStructureCheck}()},
\code{\link{.removeDisconnected}()},
\code{\link{checkFam}()}
}
\concept{check}
