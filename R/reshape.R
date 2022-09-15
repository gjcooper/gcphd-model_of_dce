#' Reshape a correlation matrix into a data.frame
#'
#' This function takes a square matrix, like a correlation matrix derived from a
#' mean covariance matrix, and reshapes it into a data.frame for plotting.
#'
#' @param x A matrix with named x and y dimensions.
#'
#' @return A df with X, Y and value columns
#'
#' @importFrom tidyr pivot_longer
#' @export
matrix_to_df <- function(x) {
  data.frame(x) %>%
    mutate(X = names(.)) %>%
    pivot_longer(-X, names_to = "Y", values_to = "value")
}
