#
build_model_name = function(model_type, test_type, lstm_type, data_type, epochs) {
  paste(model_type, lstm_type, data_type, test_type, epochs, "epochs", sep="_")
}

#normalize <- function(m){
#  (m - min(m))/(max(m)-min(m))
#}

normalize <- function(vec, min, max) {
  (vec-min) / (max-min)
}
denormalize <- function(vec,min,max) {
  vec * (max - min) + min
}

# get data into "timesteps form": single matrix, for later chop-up into X and Y parts
build_matrix <- function(tseries, overall_timesteps) {
  X <- t(sapply(1:(length(tseries) - overall_timesteps + 1), #!!!!!!!!!! +1
             function(x) tseries[x:(x + overall_timesteps - 1)]))
  cat("\nBuilt matrix with dimensions: ", dim(X))
  return(X)
}

# get data into "timesteps form": design matrix
build_X <- function(tseries, lstm_num_timesteps) {
  X <- if (lstm_num_timesteps > 1) {
    t(sapply(1:(length(tseries) - lstm_num_timesteps),
             function(x) tseries[x:(x + lstm_num_timesteps - 1)]))
  } else {
    tseries[1:length(tseries) - lstm_num_timesteps]
  }
  if (lstm_num_timesteps == 1) dim(X) <- c(length(X),1)
  cat("\nBuilt X matrix with dimensions: ", dim(X))
  return(X)
}

# get data into "timesteps form": target
build_y <- function(tseries, lstm_num_timesteps) {
  y <- sapply((lstm_num_timesteps + 1):(length(tseries)), function(x) tseries[x])
  cat("\nBuilt y vector with length: ", length(y))
  return(y)
}

# Keras LSTMs expect the input array to be shaped as (no. samples, no. time steps, no. features)
reshape_X_3d <- function(X) {
  dim(X) <- c(dim(X)[1], dim(X)[2], 1)
  cat("\nReshaped X to dimensions: ", dim(X))
  return(X)
}

# multistep forecasts, as per https://robjhyndman.com/hyndsight/rolling-forecasts/ 
# 2 variants:
# - reestimate model as new data point comes in
# - re-select complete model as new data point comes in 
# we keep the complete training set (as would be realistic)
forecast_rolling <- function(fit, n_forecast, train, test, fmode = "reestimate_only") {
  
  n <- length(test) - n_forecast + 1
  order <- arimaorder(fit)
  predictions <- matrix(0, nrow=n, ncol= n_forecast)
  lower <- matrix(0, nrow=n, ncol= n_forecast) 
  upper <- matrix(0, nrow=n, ncol= n_forecast)
  
  for(i in 1:n) {  
    x <- c(train, test[0:(i-1)])
    if(fmode == "reestimate_only") {  # re-estimate parameters, given model 
      # important: must also pass in the period because this information gets lost when converting ts to vectors in concatenation step above
      if(!is.na(order[7])) {
        refit <- Arima(x, order=order[1:3],  seasonal=list(order = order[4:6], period = order[7]))
      } else {
        refit <- Arima(x, order=order[1:3],  seasonal = order[4:6])
      }
    } else if (fmode == "recompute_model") { # re-select the whole model
      refit <- auto.arima(x)
    }
    predictions[i,] <- forecast(refit, h = n_forecast)$mean
    lower[i,] <- unclass(forecast(refit, h = n_forecast)$lower)[,2] # 95% prediction interval
    upper[i,] <- unclass(forecast(refit, h = n_forecast)$upper)[,2] # 95% prediction interval
  }
  
  list(predictions = predictions, lower = lower, upper = upper)
}

rmse <- function(target, predicted) sqrt(mean((target - predicted)^2))


