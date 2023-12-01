# Example for plots in R

# bar plot
data <- data.frame(
  continent = c("Afrika", "Antarktis", "Asien", "Australien"),
  landmass = c(30.37, 13.66, 44.58, 8.53)
)

png("~/Downloads/barplot.png", antialias = "gray")
barplot(
  height = data$landmass, names = data$continent,
  ylab = "Landmasse in Mio km2",
  cex.lab = 1.5, cex.names = 1.5
)
dev.off()

# histogram plot
data <- rnorm(1000, mean = 0, sd = 1)

png("~/Downloads/hist10.png", antialias = "gray")
hist(data,
  xlab = "Normalverteilung",
  ylab = "Häufigkeit",
  breaks = 10,
  cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5
)
dev.off()

png("~/Downloads/hist30.png", antialias = "gray")
hist(data,
  xlab = "Normalverteilung",
  ylab = "Häufigkeit",
  breaks = 30,
  cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5
)
dev.off()

png("~/Downloads/hist100.png", antialias = "gray")
hist(data,
  xlab = "Normalverteilung",
  ylab = "Häufigkeit",
  breaks = 100,
  cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5
)
dev.off()


# Boxplot of the Iris data set broken down by species
data(iris)
png("~/Downloads/boxplot.png", antialias = "gray")
boxplot(iris$Sepal.Length ~ iris$Species,
  xlab = "Species",
  ylab = "Sepal Length",
  cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.5
)
dev.off()

# Scatter plot of the Iris data set
data(iris)
png("~/Downloads/scatterplot.png", antialias = "gray")
plot(iris$Sepal.Length, iris$Sepal.Width,
  xlab = "Sepal Length",
  ylab = "Sepal Width",
  cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2
)
dev.off()

# Scatter plot of a two dimensional Gaussian with slight correlation
data <- MASS::mvrnorm(10000, mu = c(0, 0), Sigma = matrix(c(1, 0.4, 0.4, 1), nrow = 2))
png("~/Downloads/scatterplot_norm1.png", antialias = "gray")
plot(data[, 1], data[, 2],
  xlab = "X",
  ylab = "Y",
  pch = 19,
  cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2
)
dev.off()

png("~/Downloads/scatterplot_norm2.png", antialias = "gray")
plot(data[, 1], data[, 2],
  xlab = "X",
  ylab = "Y",
  pch = 19,
  col = rgb(0, 0, 1, 0.1),
  cex.lab = 1.5, cex.main = 1.5, cex.axis = 1.2
)
dev.off()