% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/data.R
\docType{data}
\name{classes}
\alias{classes}
\title{Polish school classes data set.}
\format{
An \code{mlnet} object.
}
\usage{
data(classes)
}
\description{
The Polish school classes data set \code{classes} is a subset of a larger data set which was generated as part of a Polish study on adolescent youth. The network data, obtained via a nomination processes, results in a binary, directed random graph where a directed edge from i to j indicates that student i nominated student j as a playmate. A further description of the data as well as a demonstration of an analysis with curved ERGMs can be found in Stewart, Schweinberger, Bojanowski, and Morris (2018).
}
\details{
A dataset containing network data for 9 school classes as part of a Polish educational study. 
The nodes of the network are students with nodal covariate \code{sex} and known
class membership of the students.
}
\references{
Dolata, R. (ed). (2014). 
Czy szkoła ma znaczenie? Zróżnicowanie wyników nauczania po pierwszym etapie edukacyjnym oraz jego pozaszkolne i szkolne uwarunkowania. 
Vol. 1. Warsaw: Instytut Badań Edukacyjnych.

 Dolata, R. and Rycielski, P. (2014). 
 Wprowadzenie: założenia i cele badania szkolnych uwarunkowań efektywności kształcenia SUEK.

Stewart, J., Schweinberger, M., Bojanowski, M., and M. Morris (2019).
Multilevel network data facilitate statistical inference for curved {ERGM}s with geometrically weighted terms.
Social Networks, to appear.
}
\keyword{datasets}
