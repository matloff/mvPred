source("R/utils.R") 

# Compare multiple methods using bootstrap() and return a summary table.                                                                                                                                                 
lm_compare <- function(                                                                                                                                                                                                  
  data,                                                                                                                                                                                                                  
  predictors = NULL,                                                                                                                                                                                                     
  yName,                                                                                                                                                                                                                 
  k = 5,                                                                                                                                                                                                                 
  task = c("regression", "classification"),                                                                                                                                                                              
  methods = c("CC", "AC", "TOWER", "PREFILL"),                                                                                                                                                                           
  seed = 42,                                                                                                                                                                                                             
  # PREFILL options                                                                                                                                                                                                      
  impute_method = c("mice", "amelia", "missforest", "complete"),                                                                                                                                                         
  mice_method = NULL,                                                                                                                                                                                                    
  m = 5,                                                                                                                                                                                                                 
  use_dummies = FALSE,                                                                                                                                                                                                   
  # TOWER options                                                                                                                                                                                                        
  tower_regFtnName = "lm",                                                                                                                                                                                               
  tower_opts = list(),                                                                                                                                                                                                   
  tower_scaling = NULL,                                                                                                                                                                                                  
  tower_yesYVal = NULL,                                                                                                                                                                                                  
  ...                                                                                                                                                                                                                    
) {                                                                                                                                                                                                                      
  task <- match.arg(task)                                                                                                                                                                                                
  impute_method <- match.arg(impute_method)                                                                                                                                                                              
                                                                                                                                                                                                                          
  df <- as.data.frame(data)                                                                                                                                                                                              
  if (!(yName %in% names(df))) stop(paste("Column", yName, "not found in data."))                                                                                                                                        
                                                                                                                                                                                                                          
  if (is.null(predictors)) {                                                                                                                                                                                             
    predictors <- setdiff(names(df), yName)                                                                                                                                                                              
  }                                                                                                                                                                                                                      
  missing_preds <- setdiff(predictors, names(df))                                                                                                                                                                        
  if (length(missing_preds) > 0) {                                                                                                                                                                                       
    stop(paste("Predictors not found in data:", paste(missing_preds, collapse = ", ")))                                                                                                                                  
  }                                                                                                                                                                                                                      
                                                                                                                                                                                                                          
  df <- df[, c(yName, predictors), drop = FALSE]                                                                                                                                                                         
                                                                                                                                                                                                                          
  # Drop rows with missing y for ALL methods (matches benchmark behavior)                                                                                                                                                
  df <- df[!is.na(df[[yName]]), , drop = FALSE]                                                                                                                                                                          
                                                                                                                                                                                                                          
  if (nrow(df) < k) {                                                                                                                                                                                                    
    if (task == "classification") {                                                                                                                                                                                      
      return(data.frame(                                                                                                                                                                                                 
        method = methods,                                                                                                                                                                                                
        accuracy_mean = NA_real_,                                                                                                                                                                                        
        precision_mean = NA_real_,                                                                                                                                                                                       
        recall_mean = NA_real_,                                                                                                                                                                                          
        f1_mean = NA_real_,                                                                                                                                                                                              
        auc_mean = NA_real_,                                                                                                                                                                                             
        n = nrow(df),                                                                                                                                                                                                    
        stringsAsFactors = FALSE                                                                                                                                                                                         
      ))                                                                                                                                                                                                                 
    }                                                                                                                                                                                                                    
    return(data.frame(                                                                                                                                                                                                   
      method = methods,
      MSE_mean = NA_real_,                                                                                                                                                                                               
      RMSE_mean = NA_real_,                                                                                                                                                                                              
      MAE_mean = NA_real_,                                                                                                                                                                                               
      R2_mean = NA_real_,                                                                                                                                                                                                
      n = nrow(df),                                                                                                                                                                                                      
      stringsAsFactors = FALSE                                                                                                                                                                                           
    ))                                                                                                                                                                                                                   
  }                                                                                                                                                                                                                      
                                                                                                                                                                                                                          
  dots <- list(...)                                                                                                                                                                                                      
                                                                                                                                                                                                                          
  run_one <- function(method) {                                                                                                                                                                                          
    args <- list(                                                                                                                                                                                                        
      data = df,                                                                                                                                                                                                         
      yName = yName,                                                                                                                                                                                                     
      k = k,                                                                                                                                                                                                             
      task = task,                                                                                                                                                                                                       
      method = method,                                                                                                                                                                                                   
      seed = seed                                                                                                                                                                                                        
    )                                                                                                                                                                                                                    
                                                                                                                                                                                                                          
    if (method == "PREFILL") {                                                                                                                                                                                           
      args$impute_method <- impute_method                                                                                                                                                                                
      args$mice_method <- mice_method                                                                                                                                                                                    
      args$m <- m                                                                                                                                                                                                        
      args$use_dummies <- use_dummies                                                                                                                                                                                    
    }                                                                                                                                                                                                                    
                                                                                                                                                                                                                          
    if (method == "TOWER") {                                                                                                                                                                                             
      args$tower_regFtnName <- tower_regFtnName                                                                                                                                                                          
      args$tower_opts <- tower_opts                                                                                                                                                                                      
      args$tower_scaling <- tower_scaling                                                                                                                                                                                
      args$tower_yesYVal <- tower_yesYVal                                                                                                                                                                                
    }                                                                                                                                                                                                                    
                                                                                                                                                                                                                          
    if (length(dots) > 0) {                                                                                                                                                                                              
      args <- c(args, dots)                                                                                                                                                                                              
    }                                                                                                                                                                                                                    
                                                                                                                                                                                                                          
    tryCatch(                                                                                                                                                                                                            
      do.call(bootstrap, args),                                                                                                                                                                                          
      error = function(e) {
        warning(sprintf("Method %s failed: %s", method, e$message))                                                                                                                                                      
        NULL                                                                                                                                                                                                             
      }                                                                                                                                                                                                                  
    )                                                                                                                                                                                                                    
  }                                                                                                                                                                                                                      
                                                                                                                                                                                                                          
  res_list <- lapply(methods, run_one)                                                                                                                                                                                   
                                                                                                                                                                                                                          
  if (task == "classification") {                                                                                                                                                                                        
    out <- data.frame(                                                                                                                                                                                                   
      method = methods,                                                                                                                                                                                                  
      accuracy_mean = vapply(res_list, function(r) if (is.null(r)) NA_real_ else r$accuracy_mean, numeric(1)),                                                                                                           
      precision_mean = vapply(res_list, function(r) if (is.null(r)) NA_real_ else r$precision_mean, numeric(1)),                                                                                                         
      auc_mean = vapply(res_list, function(r) if (is.null(r)) NA_real_ else r$auc_mean, numeric(1)),
      n = nrow(df),
      stringsAsFactors = FALSE
    )
  } else {
    out <- data.frame(
      method = methods,
      MSE_mean = vapply(res_list, function(r) if (is.null(r)) NA_real_ else r$MSE_mean, numeric(1)),
      RMSE_mean = vapply(res_list, function(r) if (is.null(r)) NA_real_ else r$RMSE_mean, numeric(1)),
      MAE_mean = vapply(res_list, function(r) if (is.null(r)) NA_real_ else r$MAE_mean, numeric(1)),
      R2_mean = vapply(res_list, function(r) if (is.null(r)) NA_real_ else r$R2_mean, numeric(1)),
      n = nrow(df),
      stringsAsFactors = FALSE
    )
  }

  out
}