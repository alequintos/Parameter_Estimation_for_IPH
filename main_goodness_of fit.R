## WE NEED TO RUN FIRST main_estimation_beta_lambda.R

## After estimation, we simulate data for goodness of fit

# Load necessary library for matrix exponentiation
# install.packages("expm")
library(expm)

source("fn_simulate_PH_path.R")
source("fns_auxiliary.R")     # For g and g_inv

source("fn_discrete_path.R")   # to discretize
source("fn_bridge.R")  # bridge method

source("fn_estimate_intensity_matrix.R") # estimation of Lambda

source("fns_estimate_beta.R") # For the estimation of beta

source("fn_simulate_absorption_homogenous.R") # to simulate absorption

# install.packages("beepr")  # package to make sounds
library(beepr)  # package to make sounds

# Set seed for reproducibility
set.seed(1)  

# Parameters of simulation -----------------------------------------------------

T <- 10000000

N <- 1000 # Number of trajectories 

# Define parameters of IPH
Lambda2 <- Lambda_hat_list[which(number_betas==2)[1]][[1]]
beta2 <- beta_hat_list[which(number_betas==2)[1]]


# No need to modify from here onwards ------------------------------------------
num_phases <- length(alpha)   # Number of transient phases (excluding the absorbing state)

# Generating MJP (Inhomogeneous case) --------------------------------------------

# Generate N sample paths
trajectory_IPH_real <- vector("list", 2)      # All the simulation is stored in trajectory_PH
trajectory_IPH_real[[1]] <- vector("list", N)
trajectory_IPH_real[[2]] <- vector("list", N)
names(trajectory_IPH_real) <- c("times", "phases")

for(i in 1:N) {
  result <- fn_simulate_PH_path(Lambda, alpha, T, num_phases)  # we generate homogeneous
  trajectory_IPH_real$times[[i]] <- result$times
  trajectory_IPH_real$phases[[i]] <- result$phases
  trajectory_IPH_real$times[[i]] <- sapply(trajectory_IPH_real$times[[i]], function(x) g_beta_0(x, beta)) # We transform them to non-homogeneous
}
rm(result)

trajectory_IPH2 <- vector("list", 2)      # All the simulation is stored in trajectory_PH
trajectory_IPH2[[1]] <- vector("list", N)
trajectory_IPH2[[2]] <- vector("list", N)
names(trajectory_IPH2) <- c("times", "phases")

for(i in 1:N) {
  result <- fn_simulate_PH_path(Lambda2, alpha, T, num_phases)  # we generate homogeneous
  trajectory_IPH2$times[[i]] <- result$times
  trajectory_IPH2$phases[[i]] <- result$phases
  trajectory_IPH2$times[[i]] <- sapply(trajectory_IPH2$times[[i]], function(x) g_beta_0(x, beta2)) # We transform them to non-homogeneous
}
rm(result)


## Sample IPH
sample_IPH_real <- sapply(trajectory_IPH_real$times, function(x) tail(x, 1))
sample_IPH2 <- sapply(trajectory_IPH2$times, function(x) tail(x, 1))

### Plot
# Theoretical PDF of IPH distribution ------------------------------------------
x_values_iph <- seq(0, max(sample_IPH_real, sample_IPH2), length.out = 100)
f_iph <- numeric(length(x_values_iph))
for (i in 1:length(x_values_iph)) {
  # Compute the theoretical PDF value
  f_iph[i] <- h_beta_0(x_values_iph[i], beta) * alpha %*% expm( g_inv_beta_0(x_values_iph[i] , beta) * Lambda) %*% (-Lambda %*% matrix(1, nrow = num_phases, ncol = 1))
}

x_values_iph2 <- seq(0, max(sample_IPH_real, sample_IPH2), length.out = 100)
f_iph2 <- numeric(length(x_values_iph2))
for (i in 1:length(x_values_iph2)) {
  # Compute the theoretical PDF value
  f_iph2[i] <-h_beta_0(x_values_iph2[i], beta2) * alpha %*% expm( g_inv_beta_0(x_values_iph2[i] , beta2) * Lambda2) %*% (-Lambda2 %*% matrix(1, nrow = num_phases, ncol = 1))
}

# pdf("Rplot2.pdf",
#     width = 570/72,
#     height = 326/72,
#     pointsize = 12)

plot(x_values_iph, f_iph, type = "l", 
     main = "Theoretical vs. Estimated Density", xlab = "Time", ylab = "Density", col = "blue", lwd = 3)  # Adjust ylim to fit both densities
lines(x_values_iph2, f_iph2, type = "p", col = "red", lwd = 1)
legend("topright", legend = c("Theoretical Density", "Estimated Density"), col = c("blue", "red"), lwd = 2, bty = "n")

# dev.off()


## Plotting the histograms -----------------------------------------------------

pdf("myplot.pdf",
    width = 570/72,
    height = 326/72,
    pointsize = 12)


# Compute combined x-axis limits
x_limits <- range(c(sample_IPH_real, sample_IPH2))

# Compute histograms without plotting (to get max density)
h1 <- hist(sample_IPH_real, breaks = "fd", plot = FALSE)
h2 <- hist(sample_IPH2, breaks = "fd", plot = FALSE)

# Get the maximum density from both
y_top <- ceiling(max(c(h1$density, h2$density)) / 0.005) * 0.005

hist(sample_IPH_real, probability = TRUE, col = rgb(0, 0, 1, 0.5), border = "black", breaks = "fd",
     main = paste("Observed vs. Simulated Absorption Time Distributions\nT =", T_inhomogeneous),
     xlab = "Time", ylab = "Density", xlim = x_limits, ylim = c(0, y_top))


# Add second histogram
hist(sample_IPH2, probability = TRUE, col = rgb(1, 0, 0, 0.5), border = "black", breaks = "fd", add = TRUE)

# Add legend
legend("topright", legend = c("Observed Data", "Simulated Data"),
       fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), bty = "n")

dev.off()


## Hypothesis Test ----------------------
ks.test(sample_IPH_real,sample_IPH2) # We cannot reject
