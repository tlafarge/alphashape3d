#' Volume computation
#'
#' This function calculates the volume of the \eqn{\alpha}-shape of a given
#' sample of points in the three-dimensional space.
#'
#' The function \code{volume_ashape3d} computes the volume of the
#' \eqn{\alpha}-shape for each value of \eqn{\alpha} in
#' \code{as3d$alpha[indexAlpha]} when \code{indexAlpha} is numeric.
#'
#' If \code{indexAlpha="all"} or \code{indexAlpha="ALL"} then the function
#' computes the volume of the \eqn{\alpha}-shape for all values of \eqn{\alpha}
#' in \code{as3d$alpha}.
#'
#' @param as3d An object of class \code{"ashape3d"} that represents the
#' \eqn{\alpha}-shape of a given sample of points in the three-dimensional
#' space, see \code{\link{ashape3d}}.
#' @param byComponents Logical, if FALSE (default) \code{volume_ashape3d}
#' computes the volume of the whole \eqn{\alpha}-shape. If TRUE,
#' \code{volume_ashape3d} computes the volume of each connected component of
#' the \eqn{\alpha}-shape separately.
#' @param indexAlpha A single value or vector with the indexes of
#' \code{as3d$alpha} that should be used for the computation, see Details.
#' @return If \code{indexAlpha} is a single value then the function returns the
#' volume of the \eqn{\alpha}-shape (if the argument \code{byComponents} is set
#' to FALSE) or a vector with the volumes of each connected component of the
#' \eqn{\alpha}-shape (if the argument \code{byComponents} is set to TRUE).
#'
#' Otherwise \code{volume_ashape3d} returns a list (each object in the list as
#' described above).
#' @seealso \code{\link{ashape3d}}, \code{\link{components_ashape3d}}
#' @examples
#'
#' C1 <- matrix(runif(6000), nc = 3)
#' C2 <- matrix(runif(6000), nc = 3) + 2
#' x <- rbind(C1, C2)
#' ashape3d.obj <- ashape3d(x, alpha = 0.75)
#' plot(ashape3d.obj, byComponents = TRUE)
#'
#' # Compute the volume of the alpha-shape
#' volume_ashape3d(ashape3d.obj)
#' # Compute the volumes of the connected components of the alpha-shape
#' volume_ashape3d(ashape3d.obj, byComponents = TRUE)
#'
#' @export volume_ashape3d
volume_ashape3d <-
function (as3d, byComponents = FALSE, indexAlpha = 1)
{
    if (class(as3d) != "ashape3d") {
        cat("Argument is not of class ashape3d.\n")
        return(invisible())
    }
    tetra = as3d$tetra
    x <- as3d$x
    if (class(indexAlpha) == "character")
        if (indexAlpha == "ALL" | indexAlpha == "all")
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
    for (ii in 1:dim(tetra)[1]) {
        x1 = x[tetra[ii, 1], ]
        x2 = x[tetra[ii, 2], ]
        x3 = x[tetra[ii, 3], ]
        x4 = x[tetra[ii, 4], ]
        tetra[ii, 5] = abs((x4[1] - x1[1]) * ((x2[2] - x1[2]) *
            (x3[3] - x1[3]) - (x2[3] - x1[3]) * (x3[2] - x1[2])) +
            (x4[2] - x1[2]) * ((x2[3] - x1[3]) * (x3[1] - x1[1]) -
                (x2[1] - x1[1]) * (x3[3] - x1[3])) + (x4[3] -
            x1[3]) * ((x2[1] - x1[1]) * (x3[2] - x1[2]) - (x2[2] -
            x1[2]) * (x3[1] - x1[1])))/6
    }
    volume = NULL
    if (byComponents) {
        components = components_ashape3d(as3d, indexAlpha = indexAlpha)
        indexComponents = 0
        if (length(indexAlpha) == 1)
            components = list(components)
        for (iAlpha in indexAlpha) {
            indexComponents = indexComponents + 1
            volumeComponents = rep(0, max(components[[indexComponents]]))
            for (ii in 1:dim(tetra)[1]) {
                if (tetra[ii, 5 + iAlpha] == 1)
                  volumeComponents[components[[indexComponents]][tetra[ii,
                    1]]] = volumeComponents[components[[indexComponents]][tetra[ii,
                    1]]] + tetra[ii, 5]
            }
            if (length(indexAlpha) > 1) {
                volume = c(volume, list(volumeComponents))
            }
            else {
                volume = volumeComponents
            }
        }
    }
    else {
        for (iAlpha in indexAlpha) volume = c(volume, sum(tetra[tetra[,
            5 + iAlpha] == 1, 5]))
    }
    return(volume)
}
