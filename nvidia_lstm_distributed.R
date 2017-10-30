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

train <- window(mu, end = as.Date("2017-02-28")) %>% as.ts() %>% as.vector()
test <- window(mu, start = as.Date("2017-03-01")) %>% as.ts() %>% as.vector()

num_train <- length(train)
num_test <- length(test)
num_all <- num_train + num_test

model_exists <- FALSE

lstm_num_predictions <- 6
lstm_num_timesteps <- 6
batch_size <- 1
epochs <- 500
lstm_units <- 32
lstm_type <- "stateless"
data_type <- "data_diffed_scaled"
test_type <- "MU"
model_type <- "model_lstm_time_distributed"

model_name <- build_model_name(model_type, test_type, lstm_type, data_type, epochs)

cat("\n####################################################################################")
cat("\nRunning model: ", model_name)
cat("\n####################################################################################")

train_diff <- diff(train)
test_diff <- diff(test)

# normalize
minval <- min(train_diff)
maxval <- max(train_diff)

train_diff <- normalize(train_diff, minval, maxval)
test_diff <- normalize(test_diff, minval, maxval)

matrix_train <- build_matrix(train_diff, lstm_num_timesteps + lstm_num_predictions) 
matrix_test <- build_matrix(test_diff, lstm_num_timesteps + lstm_num_predictions) 

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

df <- cbind(test, rbind(matrix(rep(NA, 12), nrow = 12, ncol = 6), pred_test_undiff) ) %>% t()
df

# df <- data_frame(time_id = 1:149,
#                  test = test)
# for(i in seq_len(nrow(pred_test))) {
#   varname <- paste0("pred_test", i)
#   df <- mutate(df, !!varname := c(rep(NA, lstm_num_timesteps+1),
#                                   rep(NA, i-1),
#                                   pred_test_undiff[i, ],
#                                   rep(NA, num_test - lstm_num_predictions - lstm_num_timesteps-i)))
# }
# 
# df <- df %>% filter(time_id %% 10== 0)
# 
# df <- df %>% gather(key = 'type', value = 'value', test:pred_test37) %>% arrange(time_id)
# 
# ggplot(df, aes(x = time_id, y = value)) + geom_line(aes(color = type, linetype=type)) 

