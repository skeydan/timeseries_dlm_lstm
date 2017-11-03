library(ggplot2)
library(forecast)
library(dlm)
library(zoo)
library("quantmod")
library(dplyr)
library(tidyr)

reticulate::py_config()
reticulate::use_condaenv("r-tensorflow", required = TRUE)
require(keras)

source("functions.R")

amd <- getSymbols("AMD", from="2015-01-01", to="2017-10-01", auto.assign = FALSE)$AMD.Adjusted
autoplot(amd)
nvda <- getSymbols("NVDA", from="2015-01-01", to="2017-10-01", auto.assign = FALSE)$NVDA.Adjusted
autoplot(nvda)
mu <- getSymbols("MU", from="2015-01-01", to="2017-10-01", auto.assign = FALSE)$MU.Adjusted
autoplot(mu)
ibm <- getSymbols("IBM", from="2015-01-01", to="2017-10-01", auto.assign = FALSE)$IBM.Adjusted
autoplot(ibm)

df <- data.frame(index(amd), amd, nvda, mu, ibm)
(tss <- read.zoo(df))
autoplot(tss) + facet_grid(Series ~ ., scales = "free_y")

###

train <- window(mu, end = as.Date("2017-02-28")) %>% as.ts() %>% as.vector()
test <- window(mu, start = as.Date("2017-03-01")) %>% as.ts() %>% as.vector()

model_exists <- TRUE

lstm_num_timesteps <- 30
batch_size <- 1
epochs <- 50
lstm_units <- 32
model_type <- "model_lstm_simple"
lstm_type <- "stateless"
data_type <- "scaled"
test_type <- "MU"

model_name <- build_model_name(model_type, test_type, lstm_type, data_type, epochs)

cat("\n####################################################################################")
cat("\nRunning model: ", model_name)
cat("\n####################################################################################")

# normalize
minval <- min(train)
maxval <- max(train)

train <- normalize(train, minval, maxval)
test <- normalize(test, minval, maxval)

X_train <- build_X(train, lstm_num_timesteps) 
y_train <- build_y(train, lstm_num_timesteps) 

X_test <- build_X(test, lstm_num_timesteps) 
y_test <- build_y(test, lstm_num_timesteps) 

# Keras LSTMs expect the input array to be shaped as (no. samples, no. time steps, no. features)
X_train <- reshape_X_3d(X_train)
X_test <- reshape_X_3d(X_test)

num_samples <- dim(X_train)[1]
num_steps <- dim(X_train)[2]
num_features <- dim(X_train)[3]

# model
if (!model_exists) {
  set.seed(22222)
  model <- keras_model_sequential() 
  model %>% 
    layer_lstm(units = lstm_units, input_shape = c(num_steps, num_features)) %>% 
    layer_dense(units = 1) %>% 
    compile(
      loss = 'mean_squared_error',
      optimizer = 'adam'
    )
  
  model %>% summary()
  
  hist <- model %>% fit( 
    X_train, y_train, batch_size = batch_size, epochs = epochs,
    validation_data = list(X_test, y_test),
    callbacks = callback_early_stopping(patience=2)
  )
  model %>% save_model_hdf5(filepath = paste0(model_name, ".h5"))
} else {
  model <- load_model_hdf5(filepath = paste0(model_name, ".h5"))
}

pred_train <- model %>% predict(X_train, batch_size = 1)
pred_test <- model %>% predict(X_test, batch_size = 1)

pred_train <- denormalize(pred_train, minval, maxval)
pred_test <- denormalize(pred_test, minval, maxval)

test_rmse <- rmse(tail(test,length(test) - lstm_num_timesteps), pred_test)

df <- data_frame(
  time_id = 1:692,
  train_ = c(denormalize(train, minval, maxval), rep(NA, length(test))),
  test_ = c(rep(NA, length(train)), denormalize(test, minval, maxval)),
  pred_train = c(rep(NA, lstm_num_timesteps), pred_train, rep(NA, length(test))),
  pred_test = c(rep(NA, length(train)), rep(NA, lstm_num_timesteps), pred_test)
)
df <- df %>% gather(key = 'type', value = 'value', train_:pred_test)
ggplot(df, aes(x = time_id, y = value)) + geom_line(aes(color = type))

