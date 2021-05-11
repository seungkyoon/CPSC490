
library(quantreg)

# Read model in given prefix
read.model <- function(prefix) {
  X_name <- sprintf('%s_X.csv', prefix)
  X <- as.matrix(read.csv(X_name))
  W1_0_name <- sprintf('%s_W1_0.csv', prefix)
  W1_0 <- as.matrix(read.csv(W1_0_name))
  W2_0_name <- sprintf('%s_W2_0.csv', prefix)
  W2_0 <- as.matrix(read.csv(W2_0_name))
  W1_name <- sprintf('%s_W1.csv', prefix)
  W1 <- as.matrix(read.csv(W1_name))
  W2_name <- sprintf('%s_W2.csv', prefix)
  W2 <- as.matrix(read.csv(W2_name))
  X_name <- sprintf('%s_X.csv', prefix)
  X <- as.matrix(read.csv(X_name))
  return(list(X=X, W1_0=W1_0, W2_0=W2_0, W1=W1, W2=W2))
}

# Function to calculate probability of model repair success
run.success <- function(model_name, trials, delta=.025, retrain=TRUE) {
  model <- read.model(model_name)
  X <- model$X
  W1_0 <- model$W1_0
  W2_0 <- model$W2_0
  if (retrain) {
    W2_0 <- 0*W2_0
  }
  W1 <- model$W1
  W2 <- model$W2
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  p <- dim(W1)[2]
  
  eps <- seq(0, 1, by=delta)
  success <- matrix(rep(0,trials*length(eps)), trials, length(eps))
  repaired <- matrix(rep(0,trials*length(eps)), trials, length(eps))
  W1.repaired <- 0*W1
  cat(sprintf("Running %d trials on model \"%s.*\", (n, d, p)=(%d, %d, %d)\n", trials, model_name, n, d, p))
  for (t in 1:trials) {
    for (i in 1:length(eps)) {
      data <- data.frame(t(X))
      not_repaired <- 0
      for (j in 1:p) {
        corrupt <- rbinom(d, 1, prob = eps[i])
        Q <- rnorm(d, mean=1)
        theta <- as.vector(W1[,j])
        data$theta <- theta + corrupt*Q - as.vector(W1_0[,j])
        tryCatch({
          fit <- rq(theta ~ 0 + ., data = data, tau=0.5, method='fn')
          v.hat <- fit$coef
          W1.repaired[,j] <- W1_0[,j] + t(X) %*% v.hat
        },
        error = function(cond) {
          W1.repaired[,j] <- W1_0[,j]
          cat(sprintf("Trial %d, %d with epsilon %f failed\n", t, j, eps[i]))
        })
        norm <- max(abs(W1[,j] - W1.repaired[,j]))
        if (norm > 1e-5) {
          not_repaired <- not_repaired + 1
          break
        }
      }
      if (not_repaired == 0) {
        corrupt <- rbinom(p, 1, prob = eps[i])
        Q <- rnorm(p, mean=1)
        eta <- as.vector(W2)
        Xtilde <- tanh(X %*% W1.repaired)
        bata <- data.frame(t(Xtilde))
        bata$eta <- eta + corrupt*Q - as.vector(W2_0)
        tryCatch({
          fit <- rq(eta ~ 0 + ., data = bata, tau=0.5, method='fn')
          u.hat <- fit$coef
          W2.repaired <- W2_0 + t(Xtilde) %*% u.hat
        },
        error = function(cond) {
          W2.repaired <- W2_0
          cat(sprintf("Trial %d with epsilon %f failed\n", t, eps[i]))
        })
        norm <- max(abs(W2 - W2.repaired))
      }
      else {
        norm <- 1
      }
      if (norm < 1e-5) {
        success[t,i] <- 1
      }
      else {
        success[t,i] <- 0
      }
      repaired[t,i] <- p-not_repaired
    }
    if (t %% (trials/20) == 0) {
      cat(sprintf("trials=%d %s\n", t, date()))
    }
  }
  return(list(success=success, repaired=repaired, epsilon=eps))
}

# Function to calculate error in repaired model matrices W and beta
run.error <- function(model_name, trials, delta=.025, retrain=TRUE) {
  model <- read.model(model_name)
  X <- model$X
  W1_0 <- model$W1_0
  W2_0 <- model$W2_0
  if (retrain) {
    W2_0 <- 0*W2_0
  }
  W1 <- model$W1
  W2 <- model$W2
  
  n <- dim(X)[1]
  d <- dim(X)[2]
  p <- dim(W1)[2]
  
  eps <-  seq(0, 1, by=delta)
  W_error <- matrix(rep(0,trials*length(eps)), trials, length(eps))
  beta_error <- matrix(rep(NA,trials*length(eps)), trials, length(eps))
  W1.repaired <- 0*W1
  cat(sprintf("Running %d trials on model \"%s.*\", (n, d, p)=(%d, %d, %d)\n",
              trials, model_name, n, d, p))
  for (t in 1:trials) {
    for (i in 1:length(eps)) {
      data <- data.frame(t(X))
      this_norm <- 0
      for (j in 1:p) {
        corrupt <- rbinom(d, 1, prob = eps[i])
        Q <- rnorm(d, mean=1)
        theta <- as.vector(W1[,j])
        data$theta <- theta + corrupt*Q - as.vector(W1_0[,j])
        tryCatch({
          fit <- rq(theta ~ 0 + ., data = data, tau=0.5, method='fn')
          v.hat <- fit$coef
          W1.repaired[,j] <- W1_0[,j] + t(X) %*% v.hat
        },
        error = function(cond) {
          W1.repaired[,j] <- W1_0[,j]
          cat(sprintf("Trial %d, %d with epsilon %f failed\n", t, j, eps[i]))
        })
        this_norm <- this_norm + mean((W1[,j] - W1.repaired[,j])^2)
      }
      W_error[t,i] <- this_norm/p
      corrupt <- rbinom(p, 1, prob = eps[i])
      Q <- rnorm(p, mean=1)
      eta <- as.vector(W2)
      Xtilde <- tanh(X %*% W1.repaired)
      bata <- data.frame(t(Xtilde))
      bata$eta <- eta + corrupt*Q - as.vector(W2_0)
      tryCatch({
        fit <- rq(eta ~ 0 + ., data = bata, tau=0.5, method='fn')
        u.hat <- fit$coef
        W2.repaired <- W2_0 + t(Xtilde) %*% u.hat
      },
      error = function(cond) {
        W2.repaired <- W2_0
        cat(sprintf("Trial %d with epsilon %f failed\n", t, eps[i]))
      })
      beta_error[t,i] <- mean((W2 - W2.repaired)^2)
    }
    if (t %% (trials/20) == 0) {
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


num_trials <- 20

successes <- run.success("colon_retrain_62_100_200", num_trials, 0.025)
save(successes, file = "successes_100_200.RData")
errors <- run.error("colon_retrain_62_100_200", num_trials, 0.025)
save(errors, file = "errors_100_200_fix.RData")

successes <- run.success("colon_retrain_62_150_300", num_trials, 0.025)
save(successes, file = "successes_150_300.RData")
errors <- run.error("colon_retrain_62_150_300", num_trials, 0.025)
save(errors, file = "errors_150_300_fix.RData")

successes <- run.success("colon_retrain_62_250_500", num_trials, 0.025)
save(successes, file = "successes_250_500rq.RData")
errors <- run.error("colon_retrain_62_250_500", num_trials, 0.025)
save(errors, file = "errors_250_500_fix.RData")

successes <- run.success("colon_retrain_62_500_1000", num_trials, 0.025)
save(successes, file = "successes_500_1000.RData")
errors <- run.error("colon_retrain_62_500_1000", num_trials, 0.025)
save(errors, file = "errors_500_1000_fix.RData")

successes <- run.success("colon_retrain_62_1000_2000", num_trials, 0.025)
save(successes, file = "successes_1000_2000.RData")
errors <- run.error("colon_retrain_62_1000_2000", num_trials, 0.025)
save(errors, file = "errors_1000_2000_fix.RData")
