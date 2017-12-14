#' Computation of the 3D \eqn{\alpha}-shape
#'
#' Implementation in R of the \eqn{\alpha}-shape of a finite set of points in
#' the three-dimensional space. The \eqn{\alpha}-shape generalizes the convex
#' hull and allows to recover the shape of non-convex and even non-connected
#' sets in 3D, given a random sample of points taken into it. Besides the
#' computation of the \eqn{\alpha}-shape, this package provides users with
#' functions to compute the volume of the \eqn{\alpha}-shape, identify the
#' connected components and facilitate the three-dimensional graphical
#' visualization of the estimated set.
#'
#' \tabular{ll}{ Package: \tab alphashape3d\cr Version: \tab 1.2\cr Date: \tab
#' 2016-02-08\cr License: \tab GPL-2\cr LazyLoad: \tab yes\cr }
#'
#' @name alphashape3d-package
#' @docType package
#' @author Thomas Lafarge, Beatriz Pateiro-Lopez.
#'
#' Maintainers: Beatriz Pateiro-Lopez <beatriz.pateiro@@usc.es>
#' @references Edelsbrunner, H., Mucke, E. P. (1994). Three-Dimensional Alpha
#' Shapes. \emph{ACM Transactions on Graphics}, 13(1), pp.43-72.
#' @keywords package
#' @useDynLib alphashape3d, .registration = TRUE
#' @import rgl
#' @import geometry
NULL
