
library(ggplot2)

data <- read.csv('./data/mountain-chain-1.csv')

ggplot(data, aes(x, y)) + geom_point() + ggtitle("Scatterplot")

get_h_opt_with_one_out_of_sample_cv <- function(kernel, x, y) {
    #
    # Get optimum h with one-out-of-sample cross-validation.
    # This function can be used with different kernels by simply
    # passing a function through the `kernel` parameter that
    # accepts `x`, `y`, `x_new`, and `h` as parameters, where `x`
    # and `y` does not include the current sample.
    #
    h_opt <- 0
    h_cv <- 100
    for (h in seq(0.005, 1, length = 1000)) {
        h_cv_current <- 0
        for (i in 1:length(x)) {
            h_cv_current <- (
                h_cv_current +
                (y[i] - kernel(x[i], x[-i], y[-i], h))^2
            )
        }
        if (h_cv_current < h_cv) {
            h_cv <- h_cv_current
            h_opt <- h
        }
    }
    return(h_opt)
}

#
# Naradaya-Watson Non-Parametric Regression with Gaussian Kernel
#

naradaya_watson_with_gaussian_kernel <- function(x_new, x, y, h) {
    d <- (x_new - x) / h
    k <- dnorm(d)
    v <- y * k
    return(sum(v) / sum(k))
}

h_opt <- get_h_opt_with_one_out_of_sample_cv(
    naradaya_watson_with_gaussian_kernel, data$x, data$y)

y_hat <- unlist(lapply(
    data$x, naradaya_watson_with_gaussian_kernel, data$x, data$y, h_opt))

ggplot(data, aes(x, y)) +
    geom_point() +
    geom_line(aes(y = y_hat)) +
    ggtitle(paste(
        "Scatterplot with optimum Naradaya-Watson Kernel (h_opt = ",
        round(h_opt, 5), ")", sep = ""))

#
# Orthogonal Polynomial Regressor (Degree = 2)
#

y_hat <- predict(lm(data$y ~ poly(data$x, 2)))

ggplot(data, aes(x, y)) +
    geom_point() +
    geom_line(aes(y = y_hat)) +
    ggtitle("Scatterplot with Orthogonal Polynomial Regressor (Degree 2)")

#
# Orthogonal Polynomial Regressor (Degree = 3)
#

y_hat <- predict(lm(data$y ~ poly(data$x, 3)))

ggplot(data, aes(x, y)) +
    geom_point() +
    geom_line(aes(y = y_hat)) +
    ggtitle("Scatterplot with Orthogonal Polynomial Regressor (Degree 3)")

#
# Naradaya-Watson Non-Parametric Regression with Epanechnikov Kernel
#

naradaya_watson_with_epanechnikov_kernel <- function(x_new, x, y, h) {
    d <- (x_new - x) / h
    k <- numeric(length(d))
    for (i in 1:length(d)) {
        if (abs(d[i]) <= 1) {
            k[i] <- 3/4*(1 - d[i]^2)
        } else {
            k[i] <- 0
        }
    }
    v <- y * k
    if (sum(k) > 0) {
        return(sum(v) / sum(k))
    } else {
        return(0)
    }
}

h_opt <- get_h_opt_with_one_out_of_sample_cv(
    naradaya_watson_with_epanechnikov_kernel, data$x, data$y)

y_hat <- unlist(lapply(
    data$x, naradaya_watson_with_epanechnikov_kernel, data$x, data$y, h_opt))

ggplot(data, aes(x, y)) +
    geom_point() +
    geom_line(aes(y = y_hat)) +
    ggtitle(paste(
        "Scatterplot with optimum Epanechnikov Kernel"))
