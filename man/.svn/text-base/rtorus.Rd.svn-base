\name{rtorus}
\alias{rtorus}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Generate points in the torus
}
\description{
This function generates \eqn{n} random points within the torus whose minor radius is \eqn{r}, major radius is \eqn{R} and center is \eqn{ct}.
}
\usage{
rtorus(n, r, R, ct = c(0, 0, 0), rotx = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{n}{
Number of observations.
}
  \item{r}{
Minor radius (radius of the tube).
}
  \item{R}{
Major radius (distance from the center of the tube to the center of the torus).
}
  \item{ct}{
A vector with the coordinates of the center of the torus.
}
\item{rotx}{
If not NULL, a rotation through an angle \code{rotx} (in radians) about the \eqn{x}-axis is performed.
}
}
\examples{
T1 <- rtorus(2000, 0.5, 2.5)
rgl.bbox()
rgl.points(T1, col = 4)

T2 <- rtorus(2000, 0.5, 2.5, ct = c(2, 0, 0.5), rotx = pi/2)
rgl.points(T2, col = 2)
}