library(keras)
library(tensorflow)
library(reticulate)
library(quantreg)

# Define model parameters here
num_filters <- 40
kernel_size <- c(19, 19)
pool_size <- 2

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
        # print(paste(fr, row, col))
        W_i <- ((fr + row - 2)*dim(X)[3]) + col
        W_i2 <- W_i + filter_shape[2] - 1
        W_j <- ((row - 1)*output_shape[2] + col - 1)*output_shape[3] + 1
        W_j2 <- W_j + filter_shape[4] - 1
        # print(paste(W_i, W_i2, W_j, W_j2))
        W[W_i:W_i2, W_j:W_j2] <- filter_row
      }
    }
  }
  return (list("W" = W, "output_shape" = output_shape))
}

# Read a model in given prefix
read.model <- function(prefix) {
  np <- import("numpy")
  X_ <- np$load("X.npy")
  X <- array_reshape(X_, c(dim(X_)[1], prod(dim(X_)[2:4])))
  load(sprintf("%s_init.RData", prefix))
  load(sprintf("%s_weights.RData", prefix))
  built_W <- build_W(init_weights[[1]], X_)
  W1_0 <- built_W$W
  os <- built_W$output_shape
  W2_0 <- init_weights[[2]]
  W1 <- weights[[1]]
  W2 <- weights[[2]]
  return(list(X=X, X_=X_, W1_0=W1_0, W2_0=W2_0, W1=W1, W2=W2, output_shape = os))
}

# Function to calculate probability of model repair success
run.success <- function(model_name, trials, delta=.025, retrain=TRUE) {
  model <- read.model(model_name)
  X_ <- model$X_
  X <- model$X
  W1_0 <- model$W1_0
  W2_0 <- model$W2_0
  if (retrain) {
    W2_0 <- 0*W2_0
  }
  W1 <- model$W1
  uncorrupted <- build_W(W1, X_)$W
  W2 <- model$W2
  output_shape <- model$output_shape
  
  pool_flatten <- keras_model_sequential()
  pool_flatten %>% layer_reshape(output_shape, input_shape = c(dim(W1_0)[2])) %>%
    layer_average_pooling_2d(pool_size) %>%
    layer_flatten()
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  p <- dim(W1_0)[2]
  p_2 <- dim(W2)[1]
  
  eps <- seq(0, 1, by=delta)
  success <- matrix(rep(0,trials*length(eps)), trials, length(eps))
  repaired <- matrix(rep(0,trials*length(eps)), trials, length(eps))
  W1.repaired <- 0*W1_0
  cat(sprintf("Running %d trials on CNN, (n, d, p)=(%d, %d, %d)\n", trials, n, d, p))
  for (t in 1:trials) {
    for (i in 1:length(eps)) {
      data <- data.frame(t(X))
      not_repaired <- 0
      
      total_params <- prod(dim(W1))
      corrupt <- array_reshape(rbinom(total_params, 1, prob = eps[i]), dim(W1))
      Q <- array_reshape(rnorm(total_params, mean=1), dim(W1))
      corrupted <- W1 + corrupt*Q
      
      corrupted <- build_W(corrupted, X_)$W
      for (j in 1:p) {
        theta <- as.vector(corrupted[,j])
        data$theta <- theta - as.vector(W1_0[,j])
        fit <- rq(theta ~ 0 + ., data = data, tau=0.5, method='fn')
        v.hat <- fit$coef
        W1.repaired[,j] <- W1_0[,j] + t(X) %*% v.hat
        norm <- max(abs(uncorrupted[,j] - W1.repaired[,j]))
        if (norm > 5e-2) {
          print(paste(t, i, j, norm))
          not_repaired <- not_repaired + 1
          break
        }
        if (j %% 100 == 0) {
          print(paste(t, i, j, norm))
        }
      }
      if (not_repaired == 0) {
        corrupt <- rbinom(p_2, 1, prob = eps[i])
        Q <- rnorm(p_2, mean=1)
        eta <- as.vector(W2)
        out <- X %*% W1.repaired
        Xtilde <- predict(pool_flatten, out)
        bata <- data.frame(t(Xtilde))
        bata$eta <- eta + corrupt*Q - as.vector(W2_0)
        fit <- rq(eta ~ 0 + ., data = bata, tau=0.5, method='fn')
        u.hat <- fit$coef
        W2.repaired <- W2_0 + t(Xtilde) %*% u.hat
        norm <- max(abs(W2 - W2.repaired))
        print(norm)
      }
      else {
        norm <- 1
      }
      if (norm < 5e-2) {
        success[t,i] <- 1
      }
      else {
        success[t,i] <- 0
      }
      repaired[t,i] <- p-not_repaired
    }
    if (t %% (trials/10) == 0) {
      cat(sprintf("trials=%d %s\n", t, date()))
    }
  }
  return(list(success=success, repaired=repaired, epsilon=eps))
}

# Function to calculate error in repaired model matrices W and beta
run.error <- function(model_name, trials, delta=.025, retrain=TRUE) {
  model <- read.model(model_name)
  X_ <- model$X_
  X <- model$X
  W1_0 <- model$W1_0
  W2_0 <- model$W2_0
  if (retrain) {
    W2_0 <- 0*W2_0
  }
  W1 <- model$W1
  uncorrupted <- build_W(W1, X_)$W
  W2 <- model$W2
  output_shape <- model$output_shape
  
  pool_flatten <- keras_model_sequential()
  pool_flatten %>% layer_reshape(output_shape, input_shape = c(dim(W1_0)[2])) %>%
    layer_average_pooling_2d(pool_size) %>%
    layer_flatten()
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  p <- dim(W1_0)[2]
  p_2 <- dim(W2)[1]
  
  eps <-  seq(0, 1, by=delta)
  W_error <- matrix(rep(0,trials*length(eps)), trials, length(eps))
  beta_error <- matrix(rep(NA,trials*length(eps)), trials, length(eps))
  W1.repaired <- 0*W1_0
  cat(sprintf("Running %d trials on CNN, (n, d, p)=(%d, %d, %d)\n", trials, n, d, p))
  for (t in 1:trials) {
    for (i in 1:length(eps)) {
      print(paste(t, i))
      data <- data.frame(t(X))
      this_norm <- 0
      
      total_params <- prod(dim(W1))
      corrupt <- array_reshape(rbinom(total_params, 1, prob = eps[i]), dim(W1))
      Q <- array_reshape(rnorm(total_params, mean=1), dim(W1))
      corrupted <- W1 + corrupt*Q
      
      corrupted <- build_W(corrupted, X_)$W
      for (j in 1:p) {
        theta <- as.vector(corrupted[,j])
        data$theta <- theta - as.vector(W1_0[,j])
        fit <- rq(theta ~ 0 + ., data = data, tau=0.5, method='fn')
        v.hat <- fit$coef
        W1.repaired[,j] <- W1_0[,j] + t(X) %*% v.hat
        this_norm <- this_norm + mean((uncorrupted[,j] - W1.repaired[,j])^2)
      }
      W_error[t,i] <- this_norm/p
      corrupt <- rbinom(p_2, 1, prob = eps[i])
      Q <- rnorm(p_2, mean=1)
      eta <- as.vector(W2)
      out <- X %*% W1.repaired
      Xtilde <- predict(pool_flatten, out)
      bata <- data.frame(t(Xtilde))
      bata$eta <- eta + corrupt*Q - as.vector(W2_0)
      tryCatch({
        fit <- rq(eta ~ 0 + ., data = bata, tau=0.5, method='fn')
        u.hat <- fit$coef
        W2.repaired <- W2_0 + t(Xtilde) %*% u.hat
        beta_error[t,i] <- mean((W2 - W2.repaired)^2)
      },
      error = function(cond) {
        cat(sprintf("Trial %d with epsilon %f failed", t, eps[i]))
      })
    }
    if (t %% (trials/10) == 0) {
      cat(sprintf("trials=%d %s\n", t, date()))
    }
  }
  return(list(W_error=W_error, beta_error=beta_error, epsilon=eps))
}

# Simple plot of successes (from run.success)
plot.success <- function(out, add=FALSE) {
  if (add) {
    lines(out$epsilon, colMeans(out$success))
  }
  else {
    plot(out$epsilon, colMeans(out$success), 'l',
         ylab='repair probability', xlab='epsilon')
  }
  points(out$epsilon, colMeans(out$success))
}

# Simple plot of W error (from run.error)
plot.W.error <- function(out, add=FALSE) {
  if (add) {
    lines(out$epsilon, colMeans(out$W_error, na.rm = TRUE))
  }
  else {
    plot(out$epsilon, colMeans(out$W_error, na.rm = TRUE), 'l',
         ylab='neuron error', xlab='epsilon')
  }
  points(out$epsilon, colMeans(out$W_error, na.rm = TRUE))
}

# Simple plot of beta error (from run.error)
plot.beta.error <- function(out, add=FALSE) {
  if (add) {
    lines(out$epsilon, colMeans(out$beta_error, na.rm = TRUE))
  }
  else {
    plot(out$epsilon, colMeans(out$beta_error, na.rm = TRUE), 'l',
         ylab='beta error', xlab='epsilon')
  }
  points(out$epsilon, colMeans(out$beta_error, na.rm = TRUE))
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

# Run 20 trials
for (i in 1:20) {
  successes <- run.success(sprintf("cnn_%d_%d", num_filters, kernel_size[1]), 1, 0.025)
  save(successes, file = sprintf("successes_cnn_%d_%d_%d.RData", num_filters, kernel_size[1], i))
  errors <- run.ntk(sprintf("cnn_%d_%d", num_filters, kernel_size[1]), 1, 0.025)
  save(errors, file = sprintf("errors_cnn_%d_%d_%d.RData", num_filters, kernel_size[1], i))
}

