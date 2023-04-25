#' Title
#'
#' @param Z
#' @param Y
#' @param va
#' @param L
#'
#' @return
#' @export
#'
#' @examples
calculate_x <- function(Z, Y, va, L) {
  if(!is.null(dim(Y))) Y <- apply(Y, 1, sum) # if Y is matrix
  if(missing(L)) {
    # check if mass balanced
    if(!all.equal(apply(Z, 1, sum) + Y, apply(Z, 2, sum) + va)) {
      stop("IO system is not mass balanced !!")
    }
    # calculate output
    x <- apply(Z, 1, sum) + Y
  } else {
    x <- L %*% Y
  }
  return(x)
}


#' Title
#'
#' @param Z
#' @param x
#'
#' @return
#' @export
#'
#' @examples
calculate_A <- function(Z, x) {
  # calculate A-matrix
  # A <- Z/x[col(Z)]
  A <- Rfast::eachrow(Z, x,'/')
  A[is.na(A)] <- 0
  return(A)
}
#' Title
#'
#' @param A
#'
#' @return
#' @export
#'
#' @examples
calculate_L <- function(A) {
  # calculate Leontief inverse
  L <- solve(diag(nrow(A)) - A)
  return(L)
}


#' Title
#' #todo> check
#' @param Z
#'
#' @return
#' @export
#'
#' @examples
calculate_B <- function(Z, x) {
  B <- Z / x
  B[is.na(B)] <- 0
  return(B)
}

#' Title
#' #TODO: check
#' @param Z
#'
#' @return
#' @export
#'
#' @examples
calculate_G <- function(B) {
  return(solve(diag(nrow(B)) -(B)))
}

#' Title
#'
#' @param E
#' @param x
#'
#' @return
#' @export
#'
#' @examples
calculate_S <- function(E, x) {
  # calculate Stressor matrix
  x_hat <- diag(1/x)
  x_hat[is.infinite(x_hat)] <- 0
  S <- E %*% x_hat
  return(S)
}

#' Title
#'
#' @param Z
#' @param Y
#' @param va
#' @param E
#'
#' @return
#' @export
#'
#' @examples
IO_creator <- function(Z, Y, va, E) {
  x <- calculate_x(Z, Y, va)
  A <- calculate_A(Z, x)
  S <- calculate_S(E, x)
  L <- calculate_L(A)
  return(list("A" = A, "L" = L, "S" = S))
}

#' Title
#'
#' @param S
#' @param L
#' @param Y
#' @param B
#' @param d
#' @param f
#' @param detailed
#'
#' @return
#' @export
#'
#' @examples
IO_calculator <- function(S, L, Y, B, d, f, detailed = TRUE) {
  if(missing(Y)) Y <- (B %*% d) * as.numeric(f)
  x <- as.numeric(L %*% Y)
  if(detailed) B <- S %*% diag(x)
  else B <- S %*% x
  return(B)
}

#' Title
#'
#' @param n.industries
#' @param n.emissions
#' @param n.fdcats
#' @param A
#'
#' @return
#' @export
#'
#' @examples
create_random_IOtable <- function(n.industries, n.emissions, n.fdcats, A = FALSE) {
  x0 <- list("S" = matrix(runif(n.industries*n.emissions), n.emissions, n.industries),
             "L" = matrix(runif(n.industries^2), n.industries, n.industries),
             "Y" = matrix(runif(n.industries * n.fdcats), n.industries, n.fdcats))
  if(A) x0[["A"]] <- matrix(runif(n.industries^2), n.industries, n.industries)
  return(x0)
}

#' Title
#'
#' @param A_mat
#' @param n
#'
#' @return
#' @export
#'
#' @examples
leontief_series_expansion <- function(A_mat, n) {
  list <- vector(mode = "list", length = n)
  list[[1]] <- diag(1, nrow = nrow(A_mat), ncol = ncol(A_mat))
  for(i in 2:n) {
    list[[i]] <- list[[i-1]] %*% A_mat
  }
  return(list)
}


#' Aggregates the Y matrix for specfic columns.
#' e.g. by country, final demand category
#'
#' @param Y the final demand matrix
#' @param groupings a vector with the groupings, same length as ncol(Y)
#'
#' @return
#' @export
#'
#' @examples
aggregate_Y <- function(Y, groupings) {
  if (length(groupings) != ncol(Y)) stop('groupings need to have the same length as nrow(Y)')
  grouping_levels <- unique(groupings)
  n_groups <- length(grouping_levels)
  Ynew <- matrix(0, nrow = nrow(Y), ncol = n_groups)
  colnames(Ynew) <- grouping_levels
  for (i in 1:n_groups) {
    Ynew[,i] <- rowsums(Y[, groupings == grouping_levels[i]])
  }
  return(Ynew)
}

