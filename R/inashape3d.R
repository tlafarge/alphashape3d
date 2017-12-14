#' Test of the inside of an \eqn{\alpha}-shape
#'
#' This function checks whether points are inside an \eqn{\alpha}-shape.
#'
#' The function \code{inashape3d} checks whether each point in \code{points} is
#' inside the \eqn{\alpha}-shape for each value of \eqn{\alpha} in
#' \code{as3d$alpha[indexAlpha]}.
#'
#' If \code{indexAlpha="all"} or \code{indexAlpha="ALL"} then the function
#' checks whether each point in \code{points} is inside the \eqn{\alpha}-shape
#' for all values of \eqn{\alpha} in \code{as3d$alpha}.
#'
#' @param as3d An object of class \code{"ashape3d"} that represents the
#' \eqn{\alpha}-shape of a given sample of points in the three-dimensional
#' space, see \code{\link{ashape3d}}.
#' @param indexAlpha A single value or vector with the indexes of
#' \code{as3d$alpha} that should be used for the computation, see Details.
#' @param points A 3-column matrix with the coordinates of the input points.
#' @return If \code{indexAlpha} is a single value then the function returns a
#' vector of boolean of length the number of input points. The element at
#' position \code{i} is \code{TRUE} if the point in \code{points[i,]} is inside
#' the \eqn{\alpha}-shape.
#'
#' Otherwise \code{inashape3d} returns a list of vectors of boolean values
#' (each object in the list as described above).
#' @seealso \code{\link{ashape3d}}
#' @examples
#'
#' T1 <- rtorus(2000, 0.5, 2)
#' T2 <- rtorus(2000, 0.5, 2, ct = c(2, 0, 0), rotx = pi/2)
#' x <- rbind(T1, T2)
#' ashape3d.obj <- ashape3d(x, alpha = 0.4)
#' # Random sample of points in a plane
#' points <- matrix(c(5*runif(10000) - 2.5, rep(0.01, 5000)), nc = 3)
#' in3d <- inashape3d(ashape3d.obj, points = points)
#' plot(ashape3d.obj, transparency = 0.2)
#' colors <- ifelse(in3d, "blue", "green")
#' rgl.points(points, col = colors)
#'
#' @export inashape3d
inashape3d <-
function (as3d, indexAlpha = 1, points)
{
    points = matrix(as.numeric(points), ncol = 3)
    triangles <- as3d$triang
    x <- as3d$x
    if (class(indexAlpha) == "character" & (indexAlpha == "ALL" |
        indexAlpha == "all"))
        indexAlpha = 1:length(as3d$alpha)
    if (any(indexAlpha > length(as3d$alpha)) | any(indexAlpha <=
        0)) {
        if (max(indexAlpha) > length(as3d$alpha))
            error = max(indexAlpha)
        else error = min(indexAlpha)
        stop(paste("indexAlpha out of bound : valid range = 1:",
            length(as3d$alpha), ", problematic value = ", error,
            sep = ""), call. = TRUE)
    }
    insideVect = NULL
    for (iAlpha in indexAlpha) {
        tr <- triangles[triangles[, 8 + iAlpha] == 2, c("tr1",
            "tr2", "tr3")]
        nbIntersect = intersectiontriangle(tr, x, points)
        insideVect <- c(insideVect, list(nbIntersect%%2 != 0))
    }
    if (length(indexAlpha) == 1)
        insideVect <- insideVect[[1]]
    return(insideVect)
}
