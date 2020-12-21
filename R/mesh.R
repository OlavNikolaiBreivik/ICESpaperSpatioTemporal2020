#' Create mesh
#'
#' This function creates the mesh to represent the Matern covariance strucutre
#'
#' @param  cutoff Parameter providing density of mesh
#' @param  cbound Boundary parameter
#' @return Mesh used
#'
#' @export
createMesh <- function(cutoff = NULL, cbound = NULL){
  maxEdge = c(cutoff,cutoff*4) # Longer distances between nodes outside if inner bounderary

  intPoints = constructIntPoints(confPred)$locUTM

  boundary.loc <- SpatialPoints(as.matrix(intPoints))
  boundary <- list(
    inla.nonconvex.hull(coordinates(boundary.loc), 20,resolution = 100),
    inla.nonconvex.hull(coordinates(boundary.loc), cbound))

  mesh <- inla.mesh.2d(boundary=boundary,
                       max.edge=maxEdge,
                       min.angle=c(30, 15),
                       cutoff=cutoff)
  print(paste("Mesh points:",mesh$n))
  return(list(mesh=mesh, barrier.triangles =NULL))
}
