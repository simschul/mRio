#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_A <- function(path, metadata = FALSE) {
  A <- data.table::fread(path, skip = 3, drop = c(1,2)) %>%
    as.matrix

  if (isTRUE(metadata)) {
    meta <- read_EB3_A_meta(path)
    attr(A, 'colnames') <- meta$colnames
    attr(A, 'rownames') <- meta$rownames
  }

  return(A)

}


#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_A_meta <- function(path) {
  colnames <- fread(file.path(path), nrows = 2,
                    drop = c(1,2)) %>%
    t %>%
    as.data.table %>%
    setnames(new = c('region', 'sector'))

  rownames <- fread(file.path(path), skip = 3,
                    select = c(1,2)) %>%
    setnames(new = c('region', 'sector'))

  return(list(colnames = colnames, rownames = rownames))

}



#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_Y <- function(path, metadata = FALSE) {
  Y <- data.table::fread(path, skip = 3, drop = c(1,2)) %>%
    as.matrix
  if (isTRUE(metadata)) {
    meta <- read_EB3_Y_meta(path)
    attr(Y, 'colnames') <- meta$colnames
    attr(Y, 'rownames') <- meta$rownames
  }
  return(Y)
}


#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_Y_meta <- function(path) {
  colnames <- fread(file.path(path), nrows = 2,
                    drop = c(1,2)) %>%
    t %>%
    as.data.table %>%
    setnames(new = c('region', 'category'))

  rownames <-fread(file.path(path), skip = 3,
                   select = c(1,2)) %>%
    setnames(new = c('region', 'sector'))

  return(list(colnames = colnames, rownames = rownames))

}


#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_va <- function(path, metadata = FALSE) {
  va <- fread(file.path(path),
              skip = 2, header = FALSE, drop = 1)[1:9] %>%
    as.matrix
  if (isTRUE(metadata)) {
    meta <- read_EB3_va_meta(path)
    attr(va, 'colnames') <- meta$colnames
    attr(va, 'rownames') <- meta$rownames
  }
  return(va)
}

#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_va_meta <- function(path) {
  colnames <- fread(file.path(path), nrows = 2,
                    drop = c(1,2), header = FALSE) %>%
    t %>%
    as.data.table %>%
    setnames(new = c('region', 'sector'))

  rownames <-fread(file.path(path), skip = 3,
                   select = c(1))[1:9,] %>%
    as.data.table %>%
    setnames(new = c('category'))

  return(list(colnames = colnames, rownames = rownames))

}

#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_S <- function(path, metadata = FALSE) {
  S <- fread(file.path(path),
             skip = 25, header = FALSE, drop = 1) %>%
    as.matrix
  if (isTRUE(metadata)) {
    meta <- read_EB3_S_meta(path)
    attr(S, 'colnames') <- meta$colnames
    attr(S, 'rownames') <-  meta$rownames
  }

  return(S)

}

#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_S_meta <- function(path) {
  colnames <- fread(file.path(path), nrows = 2,
                    drop = c(1), header = FALSE) %>%
    t %>%
    as.data.table %>%
    setnames(new = c('region', 'sector'))

  rownames <- fread(file.path(path), skip = 26,
                    select = c(1)) %>%
    setnames(new = c('category'))

  return(list(colnames = colnames, rownames = rownames))

}

#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_F <- function(path, metadata = FALSE) {
  Fmat <- fread(file.path(path),
             skip = 2, header = FALSE, drop = 1) %>%
    as.matrix
  if (isTRUE(metadata)) {
    meta <- read_EB3_F_meta(path)
    attr(Fmat, 'colnames') <- meta$colnames
    attr(Fmat, 'rownames') <-  meta$rownames
  }

  return(Fmat)

}

#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
read_EB3_F_meta <- function(path) {
  colnames <- fread(file.path(path), nrows = 2,
                    drop = c(1), header = FALSE) %>%
    t %>%
    as.data.table %>%
    setnames(new = c('region', 'sector'))

  rownames <- fread(file.path(path), skip = 3,
                    select = c(1)) %>%
    setnames(new = c('category'))

  return(list(colnames = colnames, rownames = rownames))

}


#' Title
#'
#' @param path
#' @param U
#' @param V
#' @param Y
#' @param va
#' @param metadata
#'
#' @return
#' @export
#'
#' @examples
parse_EB3_SUT <- function(path, U = TRUE, V = TRUE, Y = TRUE, va = TRUE,
                          metadata = FALSE,
                          path2meta = NULL) {
  SUT <- list()
  if (isTRUE(U)) SUT[['U']] <- read_EB3_A(file.path(path, "use.csv"),
                                          metadata = metadata)
  if (isTRUE(V)) SUT[['V']] <- read_EB3_A(file.path(path, "supply.csv"),
                                          metadata = metadata)
  if (isTRUE(Y)) SUT[['Y']] <- read_EB3_Y(file.path(path, "final_demand.csv"),
                                          metadata = metadata)
  if (isTRUE(va)) SUT[['va']] <- read_EB3_va(file.path(path, "value_added.csv"),
                                             metadata = metadata)

  if (isTRUE(metadata)) {
    if (is.null(path2meta)) {
      path2meta <- list.files(dirname(path), pattern = 'IOT_[[:digit:]]{4}_[ixi|pxp]',
                              full.names = TRUE) %>%
        stringr::str_subset('\\.zip$', negate = TRUE)

      if (length(path2meta) < 1) stop('please provide path2meta or set metadata to FALSE')
      path2meta <- path2meta[1]
    }
    SUT[['metadata']] <- parse_EB3_metadata(path2meta)
  }
  return(SUT)
}

#' Title
#'
#' @param path
#'
#' @return
#' @export
#'
#' @examples
parse_EB3_metadata <- function(path, indices = TRUE) {
  meta <- list()
  meta[['products']] <- fread(file.path(path, 'products.txt'))
  meta[['industries']] <- fread(file.path(path, 'industries.txt'))
  meta[['unit']] <- fread(file.path(path, 'unit.txt'))
  meta[['finaldemands']] <- fread(file.path(path, 'finaldemands.txt'))

  if (isTRUE(indices)) {
    meta[['indices_A']] <- read_EB3_A_meta(file.path(path, 'A.txt'))
    meta[['indices_Y']] <- read_EB3_Y_meta(file.path(path, 'Y.txt'))
    meta[['indices_S']] <- read_EB3_S_meta(file.path(path, 'satellite', 'F.txt'))
    meta[['indices_va']] <- read_EB3_va_meta(file.path(path, 'satellite', 'F.txt'))

  }

  return(meta)
}

