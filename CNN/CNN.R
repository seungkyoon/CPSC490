# install.packages('keras')
library(keras)
library(tensorflow)
library(reticulate)

# Define model parameters here
num_filters <- 40
kernel_size <- c(19, 19)
pool_size <- 2

# Train a convolutional network using Keras
train_model <- function(X, Y, seed = 1, lr = 0.05, epochs = 500, batch_size = 10, retrain = TRUE) {
  set_random_seed(seed)
  
  # Train model with both convolutional and dense layers trained
  model <- keras_model_sequential()
  model %>% layer_conv_2d(num_filters, kernel_size=kernel_size, input_shape=c(28, 28, 1),
                          use_bias = FALSE, trainable = TRUE,
                          kernel_initializer = initializer_random_normal(0, 0.1)) %>%
    layer_average_pooling_2d(pool_size) %>%
    layer_flatten() %>%
    layer_dense(1, activation = "tanh", use_bias = FALSE, trainable = TRUE,
                kernel_initializer = initializer_random_normal(0, 0.1)) %>% 
    compile(
      optimizer = optimizer_sgd(lr = lr),
      loss = 'mse'
    )
  init_weights <- get_weights(model)
  
  model %>% fit(X, Y, epochs=epochs, batch_size=batch_size)
  weights <- get_weights(model)
  
  if (retrain) {
    # Train model after fixing the convolutional layer and setting second layer to 0
    model_2 <- keras_model_sequential()
    model_2 %>% layer_conv_2d(num_filters, kernel_size=kernel_size, input_shape=c(28, 28, 1),
                              use_bias = FALSE, trainable = FALSE) %>%
      layer_average_pooling_2d(pool_size) %>%
      layer_flatten() %>%
      layer_dense(1, activation = "tanh", use_bias = FALSE, trainable = TRUE) %>% 
      compile(
        optimizer = optimizer_sgd(lr = lr),
        loss = 'mse'
      )
    
    weights[[2]] <- 0 * weights[[2]]
    set_weights(model_2, weights)
    
    model_2 %>% fit(X, Y, epochs=epochs, batch_size=batch_size)
    weights <- get_weights(model_2)
    
    model <- model_2
  }
  return (list("model" = model, "init_weights" = init_weights, "weights" = weights))
}

# Read in data
np <- import("numpy")
X <- np$load("X.npy")
Y <- np$load("Y.npy")

# Convert 0 labels to -1 to train with tanh
Y[Y==0] <- -1
ret <- train_model(X, Y, seed = 1, lr = 0.001, epochs = 500, batch_size = 10,
                       retrain = TRUE)
model <- ret$model


# Save model parameters for model repair
init_weights <- ret$init_weights
init_file <- sprintf("cnn_%d_%d_tanh_init.RData", num_filters, kernel_size[1])
save(init_weights, file = init_file)
weights <- ret$weights
weights_file <- sprintf("cnn_%d_%d_tanh_weights.RData", num_filters, kernel_size[1])
save(weights, file = weights_file)


# Function to build matrix form W of convolutional filters
build_W <- function(filters, X) {
  filter_shape <- dim(filters)
  
  model <- keras_model_sequential()
  model %>% layer_conv_2d(num_filters, kernel_size=kernel_size, input_shape=dim(X)[2:4],
                          use_bias = FALSE, trainable = FALSE)
  output_shape <- unlist(model$output_shape)
  
  # Build matrix that represents convolutional operation
  W <- matrix(0, prod(dim(X)[2:4]), prod(output_shape))
  for (fr in 1:filter_shape[1]) {
    filter_row <- filters[fr, , 1, ]
    for (row in 1:output_shape[1]) {
      for (col in 1:output_shape[2]) {
        W_i <- ((fr + row - 2)*dim(X)[3]) + col
        W_i2 <- W_i + filter_shape[2] - 1
        W_j <- ((row - 1)*output_shape[2] + col - 1)*output_shape[3] + 1
        W_j2 <- W_j + filter_shape[4] - 1
        W[W_i:W_i2, W_j:W_j2] <- filter_row
      }
    }
  }
  return (list("W" = W, "output_shape" = output_shape))
}


# Vectorize X to create design matrix (each pixel is a variable)
ret_2 <- build_W(weights[[1]], X)
output_shape <- ret_2$output_shape
W <- ret_2$W
X_ <- array_reshape(X, c(dim(X)[1], prod(dim(X)[2:4])))

test_out <- X_ %*% W
dim(test_out)


# Check that convolutional operation is the same as Keras' output
model_4 <- keras_model_sequential()
model_4 %>% layer_conv_2d(num_filters, kernel_size=kernel_size, input_shape=dim(X)[2:4],
                          use_bias = FALSE, trainable = FALSE) %>%
  layer_flatten()
set_weights(model_4, list(weights[[1]]))
out <- predict(model_4, X)

# Check that the largest error is still small (order of 1e-6)
max(abs(out - test_out))


# Reshape output and pass through model to ensure that prediction is the same
model_5 <- keras_model_sequential()
model_5 %>% layer_reshape(output_shape, input_shape = c(dim(test_out)[2])) %>%
  layer_average_pooling_2d(pool_size) %>%
  layer_flatten() %>%
  layer_dense(1, activation = "tanh", use_bias = FALSE)
set_weights(model_5, list(weights[[2]]))
preds <- predict(model_5, test_out)
original_preds <- predict(model, X)

# Check that the largest error is still small (order of 1e-6)
max(abs(preds - original_preds))

