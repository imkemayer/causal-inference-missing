# https://github.com/vdorie/npci/blob/master/R/util.R (2019-05-16)
"%not_in%" <- function(x, table) match(x, table, nomatch = 0) == 0

namedList <- function(...) {
  result <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(resultNames <- names(result))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  setNames(result, resultNames)
}

retargetCall <- function(call, symbol) {
  call[[1]] <- symbol
  call
}

stripCallArguments <- function(call, ...) {
  if (missing(call)) stop("call cannot be missing")
  extraArguments <- as.character(list(...))
  if (length(extraArguments) == 0) return(call)
  
  call[names(call) %not_in% extraArguments]
}

substituteTermInLanguage <- function(lang, oldTerm, newTerm)
{
  for (i in seq_along(lang)) {
    if (is.symbol(lang[[i]])) {
      if (lang[[i]] == oldTerm) lang[[i]] <- newTerm
    } else if (is.language(lang[[i]])) {
      lang[[i]] <- substituteTermInLanguage(lang[[i]], oldTerm, newTerm)
    }
  }
  lang
}

args2env <- function(parent, ...)
{
  parList <- list(...)
  substituteNames <- sapply(substitute(list(...)), deparse)[-1]
  if (is.null(resultNames <- names(parList))) resultNames <- substituteNames
  if (any(noNames <- resultNames == "")) resultNames[noNames] <- substituteNames[noNames]
  names(parList) <- resultNames
  
  list2env(parList, parent = parent)
}

# # https://github.com/vdorie/dbarts/blob/master/R/utility.R (2019-05-16)
# ## Turns data.frame w/factors into matrices of indicator variables. Differs from
# ## model.matrix as it doesn't drop columns for co-linearity even with multiple
# ## factors 
# makeModelMatrixFromDataFrame <- function(x, drop = TRUE) {
#   if (!is.data.frame(x)) stop('x is not a dataframe')
#   if (is.logical(drop) && (length(drop) != 1L || is.na(drop))) stop('when logical, drop must be TRUE or FALSE')
#   if (is.list(drop) && length(drop) != length(x)) stop('when list, drop must have length equal to x')
#   
#   result <- .Call(C_dbarts_makeModelMatrixFromDataFrame, x, drop)
#   attr(result, "term.labels") <- names(x)
#   result
# }