\name{plot.ashape3d}
\alias{plot.ashape3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot the \eqn{\alpha}-shape in 3D 
}
\description{
This function plots the \eqn{\alpha}-shape in 3D using the package \code{\link{rgl}}. 
}
\usage{
\method{plot}{ashape3d}(x, clear = TRUE, col = c(2, 2, 2), byComponents = FALSE,
 indexAlpha = 1, transparency = 1, walpha = FALSE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
An object of class \code{"ashape3d"} that represents the \eqn{\alpha}-shape of a given sample of points in the three-dimensional space, see \code{\link{ashape3d}}.
}
  \item{clear}{
Logical, specifying whether the current rgl device should be cleared.
}
  \item{col}{
A vector of length three specifying the colors of the triangles, edges and vertices composing the \eqn{\alpha}-shape, respectively.
}
  \item{byComponents}{
Logical, if TRUE the connected components of the \eqn{\alpha}-shape are represented in different colors, see \code{\link{components_ashape3d}}.
}
  \item{indexAlpha}{
A single value or vector with the indexes of \code{x$alpha} that should be used for the computation, see Details.
}
  \item{transparency}{
The coefficient of transparency, from 0 (transparent) to 1 (opaque), used to plot the \eqn{\alpha}-shape.  
}
  \item{walpha}{
Logical, if TRUE the value of \eqn{\alpha} is displayed in the rgl device.
}
  \item{\dots}{
Material properties. See \code{\link{rgl.material}} for details.
}
}
\details{
The function \code{plot.ashape3d} opens a rgl device for each value of \eqn{\alpha} in \code{x$alpha[indexAlpha]}. Device information is displayed in the console.

If \code{indexAlpha="all"} or \code{indexAlpha="ALL"} then the function represents the \eqn{\alpha}-shape for all values of \eqn{\alpha} in \code{as3d$alpha}.
}
\seealso{
\code{\link{ashape3d}}, \code{\link{components_ashape3d}}
}
\examples{
T1 <- rtorus(1000, 0.5, 2)
T2 <- rtorus(1000, 0.5, 2, ct = c(2, 0, 0), rotx = pi/2)
x <- rbind(T1, T2)
alpha <- c(0.15, 0.25, 1)
ashape3d.obj <- ashape3d(x, alpha = alpha)

# Plot the alpha-shape for all values of alpha
plot(ashape3d.obj, indexAlpha = "all")

# Plot the connected components of the alpha-shape for alpha=0.25
plot(ashape3d.obj, byComponents = TRUE, indexAlpha = 2)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{package }
