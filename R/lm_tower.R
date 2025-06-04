lm_tower <- function(data, yName, regFtnName = "lm", opts = list(), scaling = NULL, yesYVal = NULL) {
  # Check toweranNA package availability
  if (!requireNamespace("toweranNA", quietly = TRUE)) {
    stop("Package 'toweranNA' is required for lm_tower(). Please install it.")
  }
  
  # Input validation
  if (!is.data.frame(data)) stop("'data' must be a data frame")
  if (!(yName %in% names(data))) stop(paste0("Response variable '", yName, "' not found in data"))
  
  # Fit tower model
  tower_obj <- toweranNA::makeTower(
    data = data,
    yName = yName,
    regFtnName = regFtnName,
    opts = opts,
    scaling = scaling,
    yesYVal = yesYVal
  )
  
  # Create lm_tower object
  obj <- list(
    data = data,
    yName = yName,
    regFtnName = regFtnName,
    tower_obj = tower_obj
  )
  class(obj) <- "lm_tower"
  obj
}

summary.lm_tower <- function(object, ...) {
  cat("lm_tower model fitted using toweranNA\n")
  cat("Response variable:", object$yName, "\n")
  cat("Regression function:", object$regFtnName, "\n")
  
  # Calculate and print number of complete cases in training data
  complete_cases <- sum(stats::complete.cases(object$data))
  cat("Number of complete cases used:", complete_cases, "\n")
  
  invisible(object)
}

predict.lm_tower <- function(object, newdata, ...) {
  if (!inherits(object, "lm_tower")) stop("Object must be of class 'lm_tower'")
  
  # Predict using toweranNA's predict method, which handles missing data internally
  preds <- predict(object$tower_obj, newx = newdata, ...)
  preds
}


#library(toweranNA)

# Example dataset with missing values
#data(airquality)
#airquality$Ozone[1:10] <- NA
#airquality$Solar.R[5:15] <- NA

# Fit tower model to predict Ozone
#lm_tower_obj <- lm_tower(airquality, yName = "Ozone")

# Summarize the fitted model
#summary(lm_tower_obj)

# Prepare new data for prediction (some rows with missing values)
#newdata <- airquality[20:30, -which(names(airquality) == "Ozone")]
#newdata[1, "Solar.R"] <- NA  # introduce missing

# Predict using tower method
#preds <- predict(lm_tower_obj, newdata = newdata)
#print(preds)
