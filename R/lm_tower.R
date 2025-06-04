# lm_tower: Build Tower method model using toweranNA::makeTower
lm_tower <- function(data, yName, regFtnName = "lm", opts = NULL, scaling = NULL, yesYVal = NULL) {
  if (!requireNamespace("toweranNA", quietly = TRUE)) {
    stop("Please install the 'toweranNA' package to use lm_tower.")
  }
  
  tower_obj <- toweranNA::makeTower(
    data = data,
    yName = yName,
    regFtnName = regFtnName,
    opts = opts,
    scaling = scaling,
    yesYVal = yesYVal
  )
  
  class(tower_obj) <- c("lm_tower", class(tower_obj))
  tower_obj
}

# Predict method for lm_tower objects
predict.lm_tower <- function(object, newdata, k = 1, ...) {
  if (!requireNamespace("toweranNA", quietly = TRUE)) {
    stop("Please install the 'toweranNA' package to use predict.lm_tower.")
  }
  
  preds <- toweranNA::predict(object, newx = newdata, k = k, ...)
  preds
}

# summary method
summary.lm_tower <- function(object, ...) {
  cat("Tower method model (class 'lm_tower')\n")
  cat("Regression function:", object$regFtnName, "\n")
  cat("Training data dimensions:", dim(object$data), "\n")
  invisible(object)
}
