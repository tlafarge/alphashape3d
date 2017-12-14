#' Normal vectors computation
#'
#' This function calculates the normal vectors of all the triangles which
#' belong to the boundary of the \eqn{\alpha}-shape.
#'
#' The function \code{surfaceNormals} computes the normal vectors of all the
#' triangles which belong to the boundary of the \eqn{\alpha}-shape for each
#' value of \eqn{\alpha} in \code{x$alpha[indexAlpha]}. The magnitude of each
#' vector is equal to the area of its associated triangle.
#'
#' If \code{indexAlpha="all"} or \code{indexAlpha="ALL"} then the function
#' computes the normal vectors for all values of \eqn{\alpha} in
#' \code{as3d$alpha}.
#'
#' @param x An object of class \code{"ashape3d"} that represents the
#' \eqn{\alpha}-shape of a given sample of points in the three-dimensional
#' space, see \code{\link{ashape3d}}.
#' @param indexAlpha A single value or vector with the indexes of
#' \code{as3d$alpha} that should be used for the computation, see Details.
#' @param display Logical, if TRUE, \code{surfaceNormals} open a new rgl device
#' and display the related \eqn{\alpha}-shape and its normals vectors.
#' @param col Color of the normal vectors.
#' @param scale Scale parameter to control the length of the surface normals,
#' only affect display.
#' @param \dots Material properties. See \code{\link{rgl.material}} for
#' details.
#' @return If \code{indexAlpha} is a single value then the function returns an
#' object of class \code{"normals"} with the following components:
#' \item{normals}{Three-column matrix with the euclidean coordinates of the
#' normal vectors.} \item{centers}{Three-column matrix with the euclidean
#' coordinates of the centers of the triangles that form the
#' \eqn{\alpha}-shape.}
#'
#' Otherwise \code{surfaceNormals} returns a list of class
#' \code{"normals-List"} (each object in the list as described above).
#' @seealso \code{\link{ashape3d}}
#' @examples
#'
#' x <- rtorus(1000, 0.5, 1)
#' alpha <- 0.25
#' ashape3d.obj <- ashape3d(x, alpha = alpha)
#' surfaceNormals(ashape3d.obj, display = TRUE)
#'
#' @export surfaceNormals
surfaceNormals <-
function (x, indexAlpha = 1, display = FALSE, col = 3, scale = 1,
    ...)
{
    as3d <- x
    tetra <- as3d$tetra
    triangles <- as3d$triang
    edges <- as3d$edge
    vertex <- as3d$vertex
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
    normals.obj = NULL
    for (iAlpha in indexAlpha) {
        tr <- triangles[triangles[, 8 + iAlpha] == 2, c("tr1",
            "tr2", "tr3")]
        te <- tetra[tetra[, 5 + iAlpha] == 1, 1:4]
        normMat <- numeric(length(tr))
        middlePoint <- numeric(length(tr))
        retour <- .C("triangleNormals", as.integer(tr), dim(tr)[1],
            as.integer(te), dim(te)[1], as.numeric(x), dim(x)[1],
            normMat, middlePoint)
        normMat = matrix(retour[[7]], ncol = 3)
        middlePoint = matrix(retour[[8]], ncol = 3)
        if (display) {
            segment = matrix(rep(0, dim(normMat)[1] * 6), ncol = 3)
            for (ii in 1:dim(normMat)[1]) {
                segment[2 * ii - 1, ] = middlePoint[ii, ]
                segment[2 * ii, ] = middlePoint[ii, ] + scale *
                  normMat[ii, ]
            }
            rgl.open()
            plot(as3d, indexAlpha = iAlpha)
            segments3d(segment, col = col, alpha = 1, ...)
        }
        normals.obj <- c(normals.obj, list(list(normals = normMat,
            centers = middlePoint)))
        class(normals.obj[[length(normals.obj)]]) = "normals"
    }
    if (length(indexAlpha) == 1) {
        normals.obj <- normals.obj[[1]]
    }
    else {
        class(normals.obj) <- "normals-List"
    }
    invisible(normals.obj)
}
