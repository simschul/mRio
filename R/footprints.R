# 2. Sectoral Footprint Functions ----------------------------------------------
#' Title
#'
#' @param S_mat
#' @param L_mat
#' @param y_vec
#' @param index
#'
#' @return
#' @export
#'
#'
#' @examples
.calc.sector.fp.direct <- function(S_mat, L_mat, y_vec, index) {
  if(missing(index)) {
    fp <- IO_calculator(S_mat, L_mat, y_vec)
  } else {
    # only for 1 sector
    fp <- S_mat[,index] %*% (L_mat[index,] %*% y_vec)
  }
  return(fp)
}

#' Title
#'
#' @param S_mat
#' @param L_mat
#' @param x
#' @param index
#'
#' @return
#' @export
#'
#' @examples
.calc.sector.fp.indirect <- function(S_mat, L_mat, x, index) {
  diag(L_mat) <- 0 # all diagonal entries (e.g. input of cars into car industry) are already considered in the direct footprint calculations
  if(missing(index)) {
    fp <- S_mat %*% L_mat %*% diag(as.numeric(x))
  } else {
    fp <- S_mat %*% L_mat[,index] %*% x[index]
  }
  return(fp)
}


#' Title
#'
#' @param L_mat
#' @param S_mat
#' @param y_vec
#' @param index the index of the sector footprints are to be calculated. If missing results for ALL sectors are returned (higher computational expenses)
#' @param detailed shall footprints be returned split up by direct + indirect emissions?
#'
#' @return
#' @export
#'
#' @examples
calc_footprint_sector <- function(L_mat, S_mat, y_vec, index,
                                  detailed = FALSE) {
  direct <- .calc.sector.fp.direct(S_mat = S_mat, L_mat = L_mat,
                                   y_vec = y_vec, index = index)
  x <- calculate_x(Y = y_vec, L = L_mat)
  indirect <- .calc.sector.fp.indirect(S_mat = S_mat, L_mat = L_mat,
                                       x = x, index = index)
  if(detailed) {
    fp <- list("direct" = direct, "indirect" = indirect)
  } else {
    fp <- direct + indirect
  }
  return(fp)
}

#' Title
#'
#' @param n number of layers. recommendation >= 8
#' @param L_mat
#' @param A_mat
#' @param y_vec
#' @param S_mat
#' @param index see ?calc_footprint_sector
#'
#' @return
#' @export
#'
#' @examples
SPA_footprint_sector <- function(n = 8, L_mat, A_mat, y_vec, S_mat, index) {
  L_series <- leontief_series_expansion(A_mat, n)
  fp <- vector(mode = "list", length = n)
  fp[[1]] <- .calc.sector.fp.direct(index = index, S_mat = S_mat,
                                    L_mat = L_mat, y_vec = y_vec)

  if(missing(index)) {
    # total output
    x <- calculate_x(L = L_mat, Y = y_vec) %>% as.numeric
    for(i in 2:n) {
      fp[[i]] <- S_mat %*% L_series[[i]] %*% diag(x)
    }
  } else {
    # output of sector i
    x <- L_mat[index,] %*% y_vec
    for(i in 2:n) {
      fp[[i]] <- S_mat %*% L_series[[i]][,index] %*% x
    }
  }
  return(fp)
}
