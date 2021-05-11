library(ggplot2)
library(gridExtra)

df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) <- c("epsilon", "model", "probability")

all_files <- c("20_11", "20_15", "20_19_tanh", "20_19", "40_19", "max_20_19")
tnames <-c("20 11x11", "20 15x15", "20 19x19", "20 19x19 (sigmoid)",
           "40 19x19", "20 19x19 (max)")

for (i in 1:length(all_files)) {
  load(sprintf("successes_cnn_%s.Rdata", all_files[i]))
  probs <- colMeans(successes$success)
  epsilon <- successes$epsilon
  tmp <- data.frame(epsilon = epsilon, model = tnames[i], probability = probs)
  
  df <- rbind(df, tmp)
}
df$model = factor(df$model)

p1 <- ggplot(df, aes(epsilon, probability, linetype = model, col = model)) +
  geom_point(shape=1) +
  geom_line(size=0.8) +
  xlab('Epsilon') + ylab('Probability of Repair') +
  theme(legend.key.width=unit(1,"cm"))


df <- data.frame(matrix(ncol = 4, nrow = 0))
colnames(df) <- c("epsilon", "model", "W.error", "beta.error")

for (i in 1:length(all_files)) {
  load(sprintf("errors_cnn_%s.Rdata", all_files[i]))
  W.error <- colMeans(errors$W_error)
  beta.error <- colMeans(errors$beta_error)
  epsilon <- errors$epsilon
  tmp <- data.frame(epsilon = epsilon, model = tnames[i],
                    W.error = W.error, beta.error = beta.error)
  
  df <- rbind(df, tmp)
}
df$model <- factor(df$model)

p2 <- ggplot(df, aes(epsilon, W.error, linetype = model, col = model)) +
  geom_point(shape=1) +
  geom_line(size=0.8) +
  xlab('Epsilon') + ylab('Error of W') +
  theme(legend.key.width=unit(1,"cm"))

p3 <- ggplot(df, aes(epsilon, beta.error, linetype = model, col = model)) +
  geom_point(shape=1) +
  geom_line(size=0.8) +
  xlab('Epsilon') + ylab('Error of beta') +
  theme(legend.key.width=unit(1,"cm"))

grid.arrange(p1, p2, p3, nrow = 1)
