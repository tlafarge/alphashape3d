\name{ashape3d}
\alias{ashape3d}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
3D \eqn{\alpha}-shape computation
}
\description{
This function calculates the 3D \eqn{\alpha}-shape of a given sample of points in the three-dimensional space for \eqn{\alpha>0}.
}
\usage{
ashape3d(x, alpha, pert = FALSE, eps = 1e-09)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A 3-column matrix with the coordinates of the input points. Alternatively, an object of class \code{"ashape3d"} can be provided, see Details.
}
  \item{alpha}{
A single value or vector of values for \eqn{\alpha}.
}
  \item{pert}{
Logical. If the input points are not in general position and \code{pert} it set to TRUE the observations are perturbed by adding random noise, see Details. 
}
  \item{eps}{
Scaling factor used for data perturbation when the input points are not in general position, see Details.
}
}
\details{
If \code{x} is an object of class \code{"ashape3d"}, then \code{ashape3d} does not recompute the 3D Delaunay triangulation (it reduces the computational cost).

If the input points \code{x} are not in general position and \code{pert} is set to TRUE, the function adds random noise to the data. The noise is generated from a normal distribution with mean zero and standard deviation \code{eps*sd(x)}.
}
\value{
An object of class \code{"ashape3d"} with the following components (see Edelsbrunner and Mucke (1994) for notation):
  \item{tetra}{For each tetrahedron of the 3D Delaunay triangulation, the matrix \code{tetra} stores the indices of the sample points defining the tetrahedron (columns 1 to 4), a value that defines the intervals for which the tetrahedron belongs to the \eqn{\alpha}-complex (column 5) and for each \eqn{\alpha} a value (1 or 0) indicating whether the tetrahedron belongs to the \eqn{\alpha}-shape (columns 6 to last).}
  \item{triang}{For each triangle of the 3D Delaunay triangulation, the matrix \code{triang} stores the indices of the sample points defining the triangle (columns 1 to 3),  a value (1 or 0) indicating whether the triangle is on the convex hull (column 4), a value (1 or 0) indicating whether the triangle is attached or unattached (column 5), values that define the intervals for which the triangle belongs to the \eqn{\alpha}-complex (columns 6 to 8) and for each \eqn{\alpha} a value (0, 1, 2 or 3) indicating, respectively, that the triangle is not in the \eqn{\alpha}-shape or it is interior, regular or singular (columns 9 to last). As defined in Edelsbrunner and Mucke (1994), a simplex in the \eqn{\alpha}-complex is interior if it does not belong to the boundary of the \eqn{\alpha}-shape.
  A simplex in the \eqn{\alpha}-complex is regular if it is part of the boundary of the \eqn{\alpha}-shape and bounds some higher-dimensional simplex in the \eqn{\alpha}-complex. A simplex in the \eqn{\alpha}-complex is singular if it is part of the boundary of the \eqn{\alpha}-shape but does not bounds any higher-dimensional simplex in the \eqn{\alpha}-complex.}
   \item{edge}{For each edge of the 3D Delaunay triangulation, the matrix \code{edge} stores the indices of the sample points defining the edge (columns 1 and 2),  a value (1 or 0) indicating whether the edge is on the convex hull (column 3), a value (1 or 0) indicating whether the edge is attached or unattached (column 4), values that define the intervals for which the edge belongs to the \eqn{\alpha}-complex (columns 5 to 7) and for each \eqn{\alpha} a value (0, 1, 2 or 3) indicating, respectively, that the edge is not in the \eqn{\alpha}-shape or it is interior, regular or singular (columns 8 to last).}
   \item{vertex}{For each sample point, the matrix \code{vertex} stores the index of the point (column 1),  a value (1 or 0) indicating whether the point is on the convex hull (column 2), values that define the intervals for which the point belongs to the \eqn{\alpha}-complex (columns 3 and 4) and for each \eqn{\alpha} a value (1, 2 or 3) indicating, respectively, if the point is interior, regular or singular (columns 5 to last).}
   \item{x}{A 3-column matrix with the coordinates of the original sample points.}
\item{alpha}{A single value or vector of values of \eqn{\alpha}.}
\item{xpert}{A 3-column matrix with the coordinates of the perturbated sample points (only when the input points are not in general position and \code{pert} is set to TRUE).}
}
\references{
Edelsbrunner, H., Mucke, E. P. (1994). Three-Dimensional Alpha Shapes. \emph{ACM Transactions on Graphics}, 13(1), pp.43-72.
}


%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
T1 <- rtorus(1000, 0.5, 2)
T2 <- rtorus(1000, 0.5, 2, ct = c(2, 0, 0), rotx = pi/2)
x <- rbind(T1, T2)
# Value of alpha
alpha <- 0.25
# 3D alpha-shape
ashape3d.obj <- ashape3d(x, alpha = alpha)
plot(ashape3d.obj)

# For new values of alpha, we can use ashape3d.obj as input (faster)
alpha <- c(0.15, 1)
ashape3d.obj <- ashape3d(ashape3d.obj, alpha = alpha)
plot(ashape3d.obj, indexAlpha = 2:3)
}

\keyword{ package }
