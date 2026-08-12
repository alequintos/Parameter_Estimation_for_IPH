## WE NEED TO RUN FIRST main_estimation_beta_lambda.R along with the Weibull function

## Improved Plots for Paper for Matrix Weibull Distribution

actual_iterations <- length(beta_hat_list[1:(which(number_betas==2)[1])]) # length of the iterations (including initialization)

## --------- beta

pdf("Rplot_beta.pdf",
    width = 570/72,
    height = 326/72,
    pointsize = 12)

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
ymin <- 2
# ymin <- min(c(beta_hat_list[2:actual_iterations], beta))
# ymax <- max(c(beta_hat_list[2:actual_iterations], beta))
ymax <- 3

plot(beta_hat_list[2:actual_iterations],  # We do not include the first element of beta list because it corresponds to initialization
     type = "p", col = "black", pch = 20, bg = "black",
     main = expression(bold("Value of " * hat(beta) * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(beta)),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,  # size of axis labels
     cex.main = 1.5,  # size of main title
     xaxt = "n",
     yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
axis(side = 2, at = seq(ymin, ymax, by = 0.2))

abline(h = beta, col = "blue", lty = 1, lwd = 2) # actual value

dev.off()


# ----- Lambda[1,1]

pdf("Rplot_lambda11.pdf",
    width = 570/72,
    height = 326/72,
    pointsize = 12)

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
lambda_vals <- sapply(Lambda_hat_list, function(x) x[1, 1])[2:(actual_iterations)] # extract the values of Lambda (we do not include the 1st element of list because it corresponds to initialization)
ymin <- -3
# ymin <- min(c(lambda_vals, Lambda[1,1]))
ymax <- -2
# ymax <- max(c(lambda_vals, Lambda[1,1]))


plot(lambda_vals,
     type = "p", col = "black", pch = 20, bg = "black",
     main = expression(bold("Value of " * hat(lambda)[11] * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(lambda)[11]),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,   
     cex.main = 1.5,   
     xaxt = "n",    
     yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
axis(side = 2, at = seq(ymin, ymax, by = 0.2))

abline(h = Lambda[1,1], col = "blue", lty = 1, lwd = 2)
dev.off()


# ----- Lambda[1,2]

pdf("Rplot_lambda12.pdf",
    width = 570/72,
    height = 326/72,
    pointsize = 12)

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
lambda_vals <- sapply(Lambda_hat_list, function(x) x[1, 2])[2:(actual_iterations)] # extract the values of Lambda (we do not include the 1st element of list because it corresponds to initialization)
ymin <- 0.075
# ymin <- min(c(lambda_vals, Lambda[1,2]))
ymax <- 0.11
# ymax <- max(c(lambda_vals, Lambda[1,2]))


plot(lambda_vals,
     type = "p", col = "black", pch = 20, bg = "black",
     main = expression(bold("Value of " * hat(lambda)[12] * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(lambda)[12]),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,   
     cex.main = 1.5,   
     xaxt = "n",    
     yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
axis(side = 2, at = seq(ymin, ymax, by = 0.005))

abline(h = Lambda[1,2], col = "blue", lty = 1, lwd = 2)

dev.off()


## ------ Lambda[2,1]

pdf("Rplot_lambda21.pdf",
    width = 570/72,
    height = 326/72,
    pointsize = 12)

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
lambda_vals <- sapply(Lambda_hat_list, function(x) x[2, 1])[2:(actual_iterations)] # extract the values of Lambda (we do not include the 1st element of list because it corresponds to initialization)
ymin <- 0.008
# ymin <- min(c(lambda_vals, Lambda[2,1]))
ymax <- 0.014
# ymax <- max(c(lambda_vals, Lambda[2,1]))


plot(lambda_vals,
     type = "p", col = "black", pch = 20, bg = "black",
     main = expression(bold("Value of " * hat(lambda)[21] * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(lambda)[21]),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,   
     cex.main = 1.5,   
     xaxt = "n",    
     yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
axis(side = 2, at = seq(ymin, ymax, by = 0.001))

abline(h = Lambda[2,1], col = "blue", lty = 1, lwd = 2)

dev.off()


## ------ Lambda[2,2]

pdf("Rplot_lambda22.pdf",
    width = 570/72,
    height = 326/72,
    pointsize = 12)

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
lambda_vals <- sapply(Lambda_hat_list, function(x) x[2, 2])[2:(actual_iterations)] # extract the values of Lambda (we do not include the 1st element of list because it corresponds to initialization)
ymin <- -0.3
# ymin <- min(c(lambda_vals, Lambda[2,2]))
ymax <- -0.09
# ymax <- max(c(lambda_vals, Lambda[2,2]))


plot(lambda_vals,
     type = "p", col = "black", pch = 20, bg = "black",
     main = expression(bold("Value of " * hat(lambda)[22] * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(lambda)[22]),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,   
     cex.main = 1.5,   
     xaxt = "n",    
     yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
axis(side = 2, at = seq(ymin, ymax, by = 0.05))

abline(h = Lambda[2,2], col = "blue", lty = 1, lwd = 2)

dev.off()