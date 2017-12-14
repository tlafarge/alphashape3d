#' Plot the \eqn{\alpha}-shape in 3D
#'
#' This function plots the \eqn{\alpha}-shape in 3D using the package
#' \code{\link{rgl}}.
#'
#' The function \code{plot.ashape3d} opens a rgl device for each value of
#' \eqn{\alpha} in \code{x$alpha[indexAlpha]}. Device information is displayed
#' in the console.
#'
#' If \code{indexAlpha="all"} or \code{indexAlpha="ALL"} then the function
#' represents the \eqn{\alpha}-shape for all values of \eqn{\alpha} in
#' \code{as3d$alpha}.
#'
#' @param x An object of class \code{"ashape3d"} that represents the
#' \eqn{\alpha}-shape of a given sample of points in the three-dimensional
#' space, see \code{\link{ashape3d}}.
#' @param clear Logical, specifying whether the current rgl device should be
#' cleared.
#' @param col A vector of length three specifying the colors of the triangles,
#' edges and vertices composing the \eqn{\alpha}-shape, respectively.
#' @param byComponents Logical, if TRUE the connected components of the
#' \eqn{\alpha}-shape are represented in different colors, see
#' \code{\link{components_ashape3d}}.
#' @param indexAlpha A single value or vector with the indexes of
#' \code{x$alpha} that should be used for the computation, see Details.
#' @param transparency The coefficient of transparency, from 0 (transparent) to
#' 1 (opaque), used to plot the \eqn{\alpha}-shape.
#' @param walpha Logical, if TRUE the value of \eqn{\alpha} is displayed in the
#' rgl device.
#' @param triangles Logical, if TRUE triangles are plotted.
#' @param edges Logical, if TRUE edges are plotted.
#' @param vertices Logical, if TRUE vertices are plotted.
#' @param \dots Material properties. See \code{\link{rgl.material}} for
#' details.
#' @seealso \code{\link{ashape3d}}, \code{\link{components_ashape3d}}
#' @importFrom grDevices rainbow
#' @examples
#'
#' T1 <- rtorus(1000, 0.5, 2)
#' T2 <- rtorus(1000, 0.5, 2, ct = c(2, 0, 0), rotx = pi/2)
#' x <- rbind(T1, T2)
#' alpha <- c(0.15, 0.25, 1)
#' ashape3d.obj <- ashape3d(x, alpha = alpha)
#'
#' # Plot the alpha-shape for all values of alpha
#' plot(ashape3d.obj, indexAlpha = "all")
#'
#' # Plot the connected components of the alpha-shape for alpha=0.25
#' plot(ashape3d.obj, byComponents = TRUE, indexAlpha = 2)
#'
#' @importFrom graphics plot
#' @export
plot.ashape3d <-
function (x, clear = TRUE, col = c(2, 2, 2), byComponents = FALSE,
    indexAlpha = 1, transparency = 1, walpha = FALSE, triangles = TRUE,
    edges = TRUE, vertices = TRUE, ...)
{
    arg.triangles <- triangles
    arg.edges <- edges
    arg.vertices <- vertices
    as3d <- x
    triangles <- as3d$triang
    edges <- as3d$edge
    vertex <- as3d$vertex
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
    if (clear) {
        rgl.clear()
    }
    if (byComponents) {
        components = components_ashape3d(as3d, indexAlpha)
        if (length(indexAlpha) == 1)
            components = list(components)
        indexComponents = 0
        for (iAlpha in indexAlpha) {
            if (iAlpha != indexAlpha[1])
                rgl.open()
            if (walpha)
                title3d(main = paste("alpha =", as3d$alpha[iAlpha]))
            cat("Device ", rgl.cur(), " : alpha = ", as3d$alpha[iAlpha],
                "\n")
            indexComponents = indexComponents + 1
            components[[indexComponents]][components[[indexComponents]] ==
                -1] = 0
            colors = c("#000000", sample(rainbow(max(components[[indexComponents]]))))
            if (arg.triangles) {
                tr <- t(triangles[triangles[, 8 + iAlpha] ==
                  2 | triangles[, 8 + iAlpha] == 3, c("tr1",
                  "tr2", "tr3")])
                if (length(tr) != 0)
                  rgl.triangles(x[tr, 1], x[tr, 2], x[tr, 3],
                    col = colors[1 + components[[indexComponents]][tr]],
                    alpha = transparency, ...)
            }
            if (arg.edges) {
                ed <- t(edges[edges[, 7 + iAlpha] == 2 | edges[,
                  7 + iAlpha] == 3, c("ed1", "ed2")])
                if (length(ed) != 0)
                  rgl.lines(x[ed, 1], x[ed, 2], x[ed, 3], col = colors[1 +
                    components[[indexComponents]][ed]], alpha = transparency,
                    ...)
            }
            if (arg.vertices) {
                vt <- t(vertex[vertex[, 4 + iAlpha] == 2 | vertex[,
                  4 + iAlpha] == 3, "v1"])
                if (length(vt) != 0)
                  rgl.points(x[vt, 1], x[vt, 2], x[vt, 3], col = colors[1 +
                    components[[indexComponents]][vt]], alpha = transparency,
                    ...)
            }
        }
    }
    else {
        for (iAlpha in indexAlpha) {
            if (iAlpha != indexAlpha[1])
                rgl.open()
            if (walpha)
                title3d(main = paste("alpha =", as3d$alpha[iAlpha]))
            cat("Device ", rgl.cur(), " : alpha = ", as3d$alpha[iAlpha],
                "\n")
            if (arg.triangles) {
                tr <- t(triangles[triangles[, 8 + iAlpha] ==
                  2 | triangles[, 8 + iAlpha] == 3, c("tr1",
                  "tr2", "tr3")])
                if (length(tr) != 0)
                  rgl.triangles(x[tr, 1], x[tr, 2], x[tr, 3],
                    col = col[1], , alpha = transparency, ...)
            }
            if (arg.edges) {
                ed <- t(edges[edges[, 7 + iAlpha] == 2 | edges[,
                  7 + iAlpha] == 3, c("ed1", "ed2")])
                if (length(ed) != 0)
                  rgl.lines(x[ed, 1], x[ed, 2], x[ed, 3], col = col[2],
                    alpha = transparency, ...)
            }
            if (arg.vertices) {
                vt <- t(vertex[vertex[, 4 + iAlpha] == 2 | vertex[,
                  4 + iAlpha] == 3, "v1"])
                if (length(vt) != 0)
                  rgl.points(x[vt, 1], x[vt, 2], x[vt, 3], col = col[3],
                    alpha = transparency, ...)
            }
        }
    }
}
