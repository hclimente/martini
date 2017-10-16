#' Subgraph of vertices with an attribute
#' 
#' @description Takes a map file and:
#' 
#' @param net An igraph network.
#' @param attr An attribute of the vertices.
#' @param values Possible values of \code{attr}
#' @return A subgraph containing only the vertices with attribute equal to any of the values in \code{values}.
#' @importFrom igraph V induced_subgraph vertex_attr %>%
subnet <- function(net, attr, values) {
  vertices <- V(net)[vertex_attr(net, attr) %in% values]
  induced_subgraph(net, vertices)
}

#' Vertices with an attribute
#' 
#' @description Takes a map file and:
#' 
#' @param net An igraph network.
#' @param attr An attribute of the vertices.
#' @param values Possible values of \code{attr}
#' @return The vertices with attribute equal to any of the values in \code{values}.
#' @importFrom igraph V vertex_attr
subvert <- function(net, attr, values) {
  V(net)[vertex_attr(net, attr) %in% values]
}