#' Connected subsets computation
#'
#' This function calculates and clusters the different connected components of
#' the \eqn{\alpha}-shape of a given sample of points in the three-dimensional
#' space.
#'
#' The function \code{components_ashape3d} computes the connected components of
#' the \eqn{\alpha}-shape for each value of \eqn{\alpha} in
#' \code{as3d$alpha[indexAlpha]} when \code{indexAlpha} is numeric.
#'
#' If \code{indexAlpha="all"} or \code{indexAlpha="ALL"} then the function
#' computes the connected components of the \eqn{\alpha}-shape for all values
#' of \eqn{\alpha} in \code{as3d$alpha}.
#'
#' @param as3d An object of class \code{"ashape3d"} that represents the
#' \eqn{\alpha}-shape of a given sample of points in the three-dimensional
#' space, see \code{\link{ashape3d}}.
#' @param indexAlpha A single value or vector with the indexes of
#' \code{as3d$alpha} that should be used for the computation, see Details.
#' @return If \code{indexAlpha} is a single value then the function returns a
#' vector \code{v} of length equal to the sample size. For each sample point
#' \code{i}, \code{v[i]} represents the label of the connected component to
#' which the point belongs (for isolated points, \code{v[i]=-1}). The labels of
#' the connected components are ordered by size where the largest one (in
#' number of vertices) gets the smallest label which is one.
#'
#' Otherwise \code{components_ashape3d} returns a list of vectors describing
#' the connected components of the \eqn{\alpha}-shape for each selected value
#' of \eqn{\alpha}.
#' @seealso \code{\link{ashape3d}}, \code{\link{plot.ashape3d}}
#' @examples
#'
#' T1 <- rtorus(1000, 0.5, 2)
#' T2 <- rtorus(1000, 0.5, 2, ct = c(2, 0, 0), rotx = pi/2)
#' x <- rbind(T1, T2)
#' alpha <- c(0.25, 2)
#' ashape3d.obj <- ashape3d(x, alpha = alpha)
#' plot(ashape3d.obj, indexAlpha = "all")
#'
#' # Connected components of the alpha-shape for both values of alpha
#' comp <- components_ashape3d(ashape3d.obj, indexAlpha = "all")
#' class(comp)
#' # Number of components and points in each component for alpha=0.25
#' table(comp[[1]])
#' # Number of components and points in each component for alpha=2
#' table(comp[[2]])
#'
#' # Plot the connected components for alpha=0.25
#' plot(ashape3d.obj, byComponents = TRUE, indexAlpha = 1)
#'
#' @export components_ashape3d
components_ashape3d <-
function (as3d, indexAlpha = 1)
{
    if (class(as3d) != "ashape3d") {
        cat("Argument is not of class ashape3d.\n")
        return(invisible())
    }
    components = NULL
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
    for (iAlpha in indexAlpha) {
        nbPoints = dim(as3d$x)[1]
        edgeSelectResult = .C("edgeselect", n = as.integer(nbPoints),
            x = as.double(as3d$x), ed = as.integer(as3d$edge[,
                1:2]), nb = as.integer(dim(as3d$edge)[1]), alpha = as.double(as3d$alpha[iAlpha]),
            nbfinal = as.integer(0), PACKAGE = "alphashape3d")
        edges = as3d$edge[, 1:2][edgeSelectResult$ed[1:edgeSelectResult$nbfinal],
            ]
        adjacency = rep(list(NULL), nbPoints)
        if (dim(edges)[1] > 0)
            for (ij in 1:dim(edges)[1]) {
                adjacency[[edges[ij, 1]]] = append(adjacency[[edges[ij,
                  1]]], edges[ij, 2])
                adjacency[[edges[ij, 2]]] = append(adjacency[[edges[ij,
                  2]]], edges[ij, 1])
            }
        visited = rep(0, nbPoints)
        for (ij in 1:nbPoints) {
            if (is.null(adjacency[[ij]]))
                visited[ij] = -1
        }
        nbPart = 0
        while (!all(as.logical(visited))) {
            nbPart = nbPart + 1
            toVisit = which(visited == 0)[1]
            while (!length(toVisit) == 0) {
                cour = toVisit[1]
                neighboor = adjacency[[cour]]
                visited[cour] = nbPart
                for (ii in neighboor) {
                  if (visited[ii] == 0) {
                    toVisit = c(toVisit, ii)
                    visited[ii] = nbPart
                  }
                }
                toVisit = toVisit[-1]
            }
        }
        label = matrix(as.numeric(as.matrix(as.data.frame(table(visited)))),
            ncol = 2)
        if (any(label[, 1] == -1))
            label = label[-which(label[, 1] == -1), , drop = FALSE]
        if (dim(label)[1] > 1) {
            label = label[order(label[, 2], decreasing = TRUE),
                ]
            visited[visited == -1] = NA
            for (ii in 1:dim(label)[1]) visited[visited == label[ii,
                1]] = -ii
            visited[is.na(visited)] = 1
            visited = -visited
        }
        if (length(indexAlpha) > 1) {
            components = c(components, list(visited))
        }
        else {
            components = visited
        }
    }
    return(components)
}
