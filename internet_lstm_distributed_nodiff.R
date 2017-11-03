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

traffic_df <- read_csv("internet-traffic-data-in-bits-fr.csv", col_names = c("hour", "bits"), skip = 1)
ggplot(traffic_df, aes(x = hour, y = bits)) + geom_line() + ggtitle("Internet traffic")

internet_train <- traffic_df$bits[1:800]
internet_test <- traffic_df$bits[801:nrow(traffic_df)]

model_exists <- FALSE

lstm_num_predictions <- 168
lstm_num_timesteps <- 168
batch_size <- 1
epochs <- 500
lstm_units <- 32
lstm_type <- "stateless"
data_type <- "data_scaled"
test_type <- "INTERNET"
model_type <- "model_lstm_time_distributed"

model_name <- build_model_name(model_type, test_type, lstm_type, data_type, epochs)

cat("\n####################################################################################")
cat("\nRunning model: ", model_name)
cat("\n####################################################################################")

# normalize
minval <- min(train)
maxval <- max(train)

train <- normalize(train, minval, maxval)
test <- normalize(test, minval, maxval)

matrix_train <- build_matrix(train, lstm_num_timesteps + lstm_num_predictions) 
matrix_test <- build_matrix(test, lstm_num_timesteps + lstm_num_predictions) 

X_train <- matrix_train[ ,1:lstm_num_timesteps]
y_train <- matrix_train[ ,(lstm_num_timesteps + 1):(lstm_num_timesteps + lstm_num_predictions)]

X_test <- matrix_test[ ,1:lstm_num_timesteps]
y_test <- matrix_test[ ,(lstm_num_timesteps + 1):(lstm_num_timesteps + lstm_num_predictions)]

# Keras LSTMs expect the input array to be shaped as (no. samples, no. time steps, no. features)
X_train <- reshape_X_3d(X_train)
X_test <- reshape_X_3d(X_test)

y_train <- reshape_X_3d(y_train)
y_test <- reshape_X_3d(y_test)

num_samples <- dim(X_train)[1]
num_steps <- dim(X_train)[2]
num_features <- dim(X_train)[3]

# model
if (!model_exists) {
  set.seed(22222)
  model <- keras_model_sequential() 
  model %>% 
    layer_lstm(units = lstm_units, input_shape = c(num_steps, num_features),
               return_sequences = TRUE) %>% 
    time_distributed(layer_dense(units = 1)) %>% 
    compile(
      loss = 'mean_squared_error',
      optimizer = 'adam'
    )
  
  model %>% summary()
  
  model %>% fit( 
    X_train, y_train, batch_size = batch_size, epochs = epochs,
    validation_data = list(X_test, y_test), callbacks = callback_early_stopping(patience=2)
  )
  model %>% save_model_hdf5(filepath = paste0(model_name, ".h5"))
} else {
  model <- load_model_hdf5(filepath = paste0(model_name, ".h5"))
}

pred_train <- model %>% predict(X_train, batch_size = 1)
pred_test <- model %>% predict(X_test, batch_size = 1)

pred_train <- denormalize(pred_train, minval, maxval)
pred_test <- denormalize(pred_test, minval, maxval)

# undiff
train_add <- train[(lstm_num_timesteps+1):(length(train)-1)]
train_add_matrix <- build_matrix(train_add, lstm_num_predictions)
pred_train_undiff <- train_add_matrix + pred_train[ , , 1]

test_add <- test[(lstm_num_timesteps+1):(length(test)-1)]
test_add_matrix <- build_matrix(test_add, lstm_num_predictions)
pred_test_undiff <- test_add_matrix + pred_test[ , , 1]

df <- data_frame(time_id = 1:149,
                 test = test)
for(i in seq_len(nrow(pred_test))) {
  varname <- paste0("pred_test", i)
  df <- mutate(df, !!varname := c(rep(NA, lstm_num_timesteps+1),
                                  rep(NA, i-1),
                                  pred_test_undiff[i, ],
                                  rep(NA, num_test - lstm_num_predictions - lstm_num_timesteps-i)))
}

ggplot(df, aes(x = time_id, y =test)) + geom_line() +
  coord_cartesian(xlim = c(45,100)) +
 geom_line(aes(y = pred_test45), color = "cyan") + 
  geom_line(aes(y = pred_test55), color = "red") + 
  geom_line(aes(y = pred_test65), color = "green") + 
  geom_line(aes(y = pred_test75), color = "violet") + 
  geom_line(aes(y = pred_test85), color = "blue") 
  

df <- df %>% gather(key = 'type', value = 'value', test:pred_test37) %>% arrange(time_id)
ggplot(df, aes(x = time_id, y = value)) + geom_line(aes(color = type, linetype=type))

