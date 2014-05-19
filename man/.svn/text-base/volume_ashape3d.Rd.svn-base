\name{volume_ashape3d}
\alias{volume_ashape3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Volume computation
}
\description{
This function calculates the volume of the \eqn{\alpha}-shape of a given sample of points in the three-dimensional space. 
}
\usage{
volume_ashape3d(as3d, byComponents = FALSE, indexAlpha = 1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{as3d}{
 An object of class \code{"ashape3d"} that represents the \eqn{\alpha}-shape of a given sample of points in the three-dimensional space, see \code{\link{ashape3d}}.
}
  \item{byComponents}{
Logical, if FALSE (default) \code{volume_ashape3d} computes the volume of the whole \eqn{\alpha}-shape. If TRUE, \code{volume_ashape3d} computes the volume of each connected component of the \eqn{\alpha}-shape separately.
}
  \item{indexAlpha}{
A single value or vector with the indexes of \code{as3d$alpha} that should be used for the computation, see Details.
}
}
\details{
The function \code{volume_ashape3d} computes the volume of the \eqn{\alpha}-shape for each value of \eqn{\alpha} in \code{as3d$alpha[indexAlpha]} when \code{indexAlpha} is numeric.

If \code{indexAlpha="all"} or \code{indexAlpha="ALL"} then the function computes the volume of the \eqn{\alpha}-shape for all values of \eqn{\alpha} in \code{as3d$alpha}.
}
\value{
If \code{indexAlpha} is a single value then the function returns the volume of the \eqn{\alpha}-shape (if the argument \code{byComponents} is set to FALSE) or a vector with the volumes of each connected component of the \eqn{\alpha}-shape (if the argument \code{byComponents} is set to TRUE). 

Otherwise \code{volume_ashape3d} returns a list (each object in the list as described above).
}
\seealso{
\code{\link{ashape3d}}, \code{\link{components_ashape3d}}
}
\examples{
C1 <- matrix(runif(6000), nc = 3)
C2 <- matrix(runif(6000), nc = 3) + 2
x <- rbind(C1, C2)
ashape3d.obj <- ashape3d(x, alpha = 0.75)
plot(ashape3d.obj, byComponents = TRUE)

# Compute the volume of the alpha-shape
volume_ashape3d(ashape3d.obj)
# Compute the volumes of the connected components of the alpha-shape
volume_ashape3d(ashape3d.obj, byComponents = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ package}
