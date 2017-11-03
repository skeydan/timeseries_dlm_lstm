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
all <- mu %>% as.ts() %>% as.vector()

(train_fit <- auto.arima(train))
(test_fit <- auto.arima(test))
(all_fit <- auto.arima(all))

num_train <- length(train)
num_test <- length(test)
num_all <- num_train + num_test

model_exists <- TRUE

lstm_num_predictions <- 30
lstm_num_timesteps <- 30
#lstm_num_predictions <- 6
#lstm_num_timesteps <- 6
#lstm_num_predictions <- 12
#lstm_num_timesteps <- 12

batch_size <- 1
epochs <- 500
lstm_units <- 32
lstm_type <- "stateless"
data_type <- "scaled"
test_type <- "MU"
model_type <- "model_lstm_time_distributed"

build_model_name = function(model_type, test_type, lstm_type, data_type, epochs) {
  paste(model_type, lstm_type, data_type, test_type, "timesteps", lstm_num_timesteps, sep="_")
}
model_name <- build_model_name(model_type, test_type, lstm_type, data_type, lstm_num_timesteps)

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

pred_train <- pred_train[ , , 1]
pred_test <- pred_test[ , , 1]

df <- data_frame(time_id = 1:length(test),
                 test = denormalize(test, minval, maxval))
for(i in seq_len(nrow(pred_test))) {
  varname <- paste0("pred_test", i)
  df <- mutate(df, !!varname := c(rep(NA, lstm_num_timesteps),
                                  rep(NA, i-1),
                                  pred_test[i, ],
                                  rep(NA, num_test - lstm_num_predictions - lstm_num_timesteps -i +1)))
}

ggplot(df, aes(x = time_id, y =test)) + geom_line() +
  coord_cartesian(xlim = c(1,75)) +
  geom_line(aes(y = pred_test8), color = "cyan") + 
  geom_line(aes(y = pred_test18), color = "red") + 
  geom_line(aes(y = pred_test28), color = "green") + 
  geom_line(aes(y = pred_test38), color = "violet") + 
  geom_line(aes(y = pred_test48), color = "blue") + 
  geom_line(aes(y = pred_test58), color = "brown") + 
  geom_line(aes(y = pred_test68), color = "orange") 

ggplot(df, aes(x = time_id, y =test)) + geom_line() +
  coord_cartesian(xlim = c(71,150)) +
  geom_line(aes(y = pred_test78), color = "cyan") + 
  geom_line(aes(y = pred_test88), color = "red") + 
  geom_line(aes(y = pred_test98), color = "violet") + 
  geom_line(aes(y = pred_test108), color = "blue") + 
  geom_line(aes(y = pred_test118), color = "green") + 
  geom_line(aes(y = pred_test128), color = "orange") + 
  geom_line(aes(y = pred_test138), color = "pink") 

ggplot(df, aes(x = time_id, y =test)) + geom_line() +
  geom_line(aes(y = pred_test6), color = "cyan") + 
  geom_line(aes(y = pred_test26), color = "red") + 
  geom_line(aes(y = pred_test46), color = "violet") + 
  geom_line(aes(y = pred_test66), color = "blue") + 
  geom_line(aes(y = pred_test86), color = "green") +
  geom_line(aes(y = pred_test106), color = "pink") +
  geom_line(aes(y = pred_test126), color = "cyan") 

ggplot(df, aes(x = time_id, y =test)) + geom_line() +
  geom_line(aes(y = pred_test1), color = "cyan") + 
  geom_line(aes(y = pred_test30), color = "violet") + 
  geom_line(aes(y = pred_test60), color = "blue") + 
  geom_line(aes(y = pred_test90), color = "red") 

