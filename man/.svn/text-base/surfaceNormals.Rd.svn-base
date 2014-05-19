\name{surfaceNormals}
\alias{surfaceNormals}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Normal vectors computation
}
\description{
This function calculates the normal vectors of all the triangles which belong to the boundary of the \eqn{\alpha}-shape. 
}
\usage{
surfaceNormals(x, indexAlpha = 1, display = FALSE, col = 3, scale = 1, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"ashape3d"} that represents the \eqn{\alpha}-shape of a given sample of points in the three-dimensional space, see \code{\link{ashape3d}}.
}
  \item{indexAlpha}{
A single value or vector with the indexes of \code{as3d$alpha} that should be used for the computation, see Details. 
}
  \item{display}{
Logical, if TRUE, \code{surfaceNormals} open a new rgl device and display the related \eqn{\alpha}-shape and its normals vectors.

}
\item{col}{
Color of the normal vectors.
}
\item{scale}{
Scale parameter to control the length of the surface normals, only affect display.
}
\item{\dots}{
Material properties. See \code{\link{rgl.material}} for details.
}
}
\details{
The function \code{surfaceNormals} computes the normal vectors of all the triangles which belong to the boundary of the \eqn{\alpha}-shape for each value of \eqn{\alpha} in \code{x$alpha[indexAlpha]}. The magnitude of each vector is equal to the area of its associated triangle.

If \code{indexAlpha="all"} or \code{indexAlpha="ALL"} then the function computes the normal vectors for all values of \eqn{\alpha} in \code{as3d$alpha}.
}
\value{
If \code{indexAlpha} is a single value then the function returns an object of class \code{"normals"} with the following components:
\item{normals}{Three-column matrix with the euclidean coordinates of the normal vectors.}
\item{centers}{Three-column matrix with the euclidean coordinates of the centers of the triangles that form the \eqn{\alpha}-shape.}

Otherwise \code{surfaceNormals} returns a list of class \code{"normals-List"} (each object in the list as described above).
}
\seealso{
\code{\link{ashape3d}}
}
\examples{
x <- rtorus(1000, 0.5, 1)
alpha <- 0.25
ashape3d.obj <- ashape3d(x, alpha = alpha)
surfaceNormals(ashape3d.obj, display = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package}
