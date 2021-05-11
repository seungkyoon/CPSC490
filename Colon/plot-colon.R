library(ggplot2)
library(gridExtra)

df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) <- c("epsilon", "p", "probability")
for (d in c(100, 150, 250, 500, 1000)) {
  p <- 2 * d
  load(sprintf("successes_%d_%d.Rdata", d, p))
  probs <- colMeans(successes$success)
  epsilon <- successes$epsilon
  tmp <- data.frame(epsilon = epsilon, p = p, probability = probs)
  print(nrow(tmp))
  
  df <- rbind(df, tmp)
}
df$p <- factor(df$p)

p1 <- ggplot(df, aes(epsilon, probability, linetype=p, col = p)) +
  geom_point(shape=1) +
  geom_line(size=0.8) +
  xlab('Epsilon') + ylab('Probability of Repair') +
  theme(legend.key.width=unit(1,"cm"))


df <- data.frame(matrix(ncol = 3, nrow = 0))
colnames(df) <- c("epsilon", "p", "W.error")
for (d in c(100, 150, 250, 500, 1000)) {
  p <- 2 * d
  load(sprintf("errors_%d_%d.Rdata", d, p))
  W.error <- colMeans(errors$W_error)
  epsilon <- errors$epsilon
  tmp <- data.frame(epsilon = epsilon, p = p, W.error = W.error)
  print(nrow(tmp))
  
  df <- rbind(df, tmp)
}
df$p <- factor(df$p)

p2 <- ggplot(df, aes(epsilon, W.error, linetype=p, col = p)) +
  geom_point(shape=1) +
  geom_line(size=0.8) +
  xlab('Epsilon') + ylab('Error of W') +
  theme(legend.key.width=unit(1,"cm"))

grid.arrange(p1, p2, nrow = 1)
