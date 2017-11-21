#' Subgraph of vertices with an attribute
#' 
#' @description Returns a subgraph matching some condition.
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
#' @description Returns the nodes matching some condition.
#' 
#' @param net An igraph network.
#' @param attr An attribute of the vertices.
#' @param values Possible values of \code{attr}
#' @return The vertices with attribute equal to any of the values in \code{values}.
#' @importFrom igraph V vertex_attr
subvert <- function(net, attr, values) {
  V(net)[vertex_attr(net, attr) %in% values]
}

#' Check package is installed
#' 
#' @description Checks if a package is installed, launches an error if it is not.
#' 
#' @param pkg Name of the package.
#' @param fn Function calling the check.
check_installed <- function(pkg, fn = "this function") {
  if (!requireNamespace(pkg, quietly = TRUE)) {
    stop(paste(pkg, "needed for", fn, "to work. Please install it."),
         call. = FALSE)
  }
}