parse.formula <- function(fm, fun) {
  l <- as.list(attr(terms(fm), "variables"))[-1]
  l[grep(fun, l)]
}