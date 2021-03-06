\name{mapvariable}
\Rdversion{1.1}
\alias{mapvariable}
\title{
Plot Levels of a Variable in a Colour-Coded Map
}
\description{
Plot levels of a variable in a colour-coded map along with a legend.
}
\usage{
mapvariable(y, spatial.polygon, ncut=1000, nlevels=10, lower=NULL, upper=NULL, 
main=NULL, xlab=NULL, ylab=NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{variable to plot}
  \item{spatial.polygon}{an object of class SpatialPolygons (See \link[sp]{SpatialPolygons-class})}
  \item{ncut}{number of cuts in colour levels to plot}
  \item{nlevels}{number of levels to include in legend}
  \item{lower}{lower bound of levels}
  \item{upper}{upper bound of levels}
  \item{main}{an overall title for the plot}
  \item{xlab}{a title for the x axis}
  \item{ylab}{a title for the y axis}
}
%\details{}
\value{
A map colour-coded to indicate the different levels of \code{y}
}
\references{
Bivand, R. S., Pebesma E. J., and Gomez-Rubio V. (2008) \emph{Applied Spatial Data Analysis with R}.  Springer Series in Statistics.

E. J. Pebesma and R. S. Bivand. (2005) Classes and methods for spatial data in R. \emph{R News}, \bold{5}, 9--13.  
}
\author{Jon Wakefield, Nicky Best, Sebastien Haneuse, and Albert Y. Kim}
%\note{}
\examples{
data(scotland)
map <- scotland$spatial.polygon
y <- scotland$data$cases
E <- scotland$data$expected
SMR <- y/E
mapvariable(SMR,map,main="Scotland",xlab="Eastings (km)",ylab="Northings (km)")
}
\keyword{file}
