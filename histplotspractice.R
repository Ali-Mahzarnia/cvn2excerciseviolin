# install.packages("vioplot")
library("vioplot")
x <- c(6, 9, 0, 19, -1, 8, 12, 5, 3, 7,
       2, 4, 3, -8, -9, 8, 4, 12, 5, 14)
vioplot(x)
vioplot(x, horizontal = TRUE)
vioplot(x,
        col = 2,               # Color of the area
        rectCol = "red",       # Color of the rectangle
        lineCol = "white",     # Color of the line
        colMed = "green",      # Pch symbol color
        border = "black",      # Color of the border of the violin
        pchMed = 16,           # Pch symbol for the median
        plotCentre = "points") # If "line", plots a median line

stripchart(x, method = "jitter", col = "blue",
           vertical = TRUE, pch = 19, add = TRUE)
box <- boxplot(x)

x <- x[!(x %in% box$out)]

vioplot(x)


par(mfrow = c(1, 2))

vioplot(1:10)
vioplot(1:10, ylog = TRUE)

par(mfrow = c(1, 1))


set.seed(1)

# Multimodal data
n <- 10000
ii <- rbinom(n, 1, 0.5)
data <- rnorm(n, mean = 130, sd = 10) * ii +
  rnorm(n, mean = 80, sd = 5) * (1 - ii)

# Histogram
hist(data, probability = TRUE, col = "grey", axes = FALSE,
     main = "", xlab = "",  ylab = "")

# X-axis
axis(1)

# Density
lines(density(data), lwd = 2, col = "red")

# Add violin plot
par(new = TRUE)
vioplot(data, horizontal = TRUE, yaxt = "n", axes = FALSE,
        col = rgb(0, 1, 1, alpha = 0.15))

tail(chickwts) # Last rows
data <- chickwts

vioplot(data$weight ~ data$feed, col = 2:length(levels(data$feed)),
        xlab = "Feed", ylab = "Weight")

stripchart(data$weight ~ data$feed, vertical = TRUE, method = "jitter",
           pch = 19, add = TRUE, col = 3:8)
