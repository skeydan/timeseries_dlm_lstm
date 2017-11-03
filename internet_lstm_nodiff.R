library(ggplot2)
library(dplyr)
library(tidyr)
library(readr)

source("common.R")
source("functions.R")

traffic_df <- read_csv("internet-traffic-data-in-bits-fr.csv", col_names = c("hour", "bits"), skip = 1)
ggplot(traffic_df, aes(x = hour, y = bits)) + geom_line() + ggtitle("Internet traffic")

internet_train <- traffic_df$bits[1:800]
internet_test <- traffic_df$bits[801:nrow(traffic_df)]

model_exists <- FALSE

lstm_num_timesteps <- 7*24
batch_size <- 1
epochs <- 500
lstm_units <- 32
model_type <- "model_lstm_simple"
lstm_type <- "stateless"
data_type <- "data_scaled"
test_type <- "INTERNET"

model_name <- build_model_name(model_type, test_type, lstm_type, data_type, epochs)

cat("\n####################################################################################")
cat("\nRunning model: ", model_name)
cat("\n####################################################################################")

train <- internet_train[!is.na(internet_train)]
test <- internet_test[!is.na(internet_test)]

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
  
  model %>% fit( 
    X_train, y_train, batch_size = batch_size, epochs = epochs,
    validation_data = list(X_test, y_test),
    callbacks = callback_early_stopping(patience=2)train = c(internet_train, rep(NA, length(internet_test))),
    
  )
  model %>% save_model_hdf5(filepath = paste0(model_name, ".h5"))
} else {
  model <- load_model_hdf5(filepath = paste0(model_name, ".h5"))
}

pred_train <- model %>% predict(X_train, batch_size = 1)
pred_test <- model %>% predict(X_test, batch_size = 1)

pred_train <- denormalize(pred_train, minval, maxval)
pred_test <- denormalize(pred_test, minval, maxval)

test_rmse <- rmse(tail(internet_test,length(internet_test) - lstm_num_timesteps), pred_test)

df <- data_frame(
                 time_id = 1:1231,
                 train = c(internet_train, rep(NA, length(internet_test))),
                 test = c(rep(NA, length(internet_train)), internet_test),
                 pred_train_ = c(rep(NA, lstm_num_timesteps), pred_train, rep(NA, length(internet_test))),
                 pred_test_ = c(rep(NA, length(internet_train)), rep(NA, lstm_num_timesteps), pred_test)
   )
df <- df %>% gather(key = 'type', value = 'value', train:pred_test_)
ggplot(df, aes(x = time_id, y = value)) + geom_line(aes(color = type))


