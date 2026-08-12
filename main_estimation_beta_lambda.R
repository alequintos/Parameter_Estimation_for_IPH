# We estimate both beta and Lambda
# In each iteration, we estimate multiple betas, but only one Lambda

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

## Yslas paper ----------------------------------------
T_inhomogeneous <- 60
# T_inhomogeneous <- 47
# T_inhomogeneous <- 39
# T_inhomogeneous <- 36

delta <- 1

# Define parameters of IPH
alpha <- c(0.0451, 0.1303, 0.8246)
Lambda <- matrix(c(-0.1357, 0.1214, 0, 0.0130, -0.0421, 0.0288, 0.1415, 0.0184, -0.1620), nrow = length(alpha), byrow = TRUE)
beta <- 0.1019

beta_0 <- 1

# Parameters for GD
eta <- 0.00000001
error_likelihood <- 0.01

## Bladt and Peralta paper ----------------------------------------
# T_inhomogeneous <- 5
# 
# delta <- 0.1
# 
# # Define parameters of IPH
# alpha <- c(0.5, 0.5)
# Lambda <- matrix(c(-3, 0.1, 0.01, -0.1), nrow = length(alpha), byrow = TRUE)
# beta <- 3
# 
# beta_0 <- 2
# 
# # Parameters for GD
# eta <- 0.001
# error_likelihood <- 0.001


## General parameters ----------------------------------------
N <- 1000 # Number of trajectories
T_homogeneous <- g_inv_beta_0(T_inhomogeneous, beta)
total_iterations <- 100  # number of iterations in estimation algorithm 
iteration<-1

Lambda_hat_list <- list()
num_phases <- length(alpha)   # Number of transient phases (excluding the absorbing state)

# Generating MJP (Inhomogeneous case) --------------------------------------------

# Preallocate variables
trajectory_IPH <- list(
  times = vector("list", length = N),
  phases = vector("list", length = N),
  discrete_times = vector("list", N),
  discrete_phases = vector("list", N)
)

# Simulate sample paths
for(i in 1:N) {
  result <- fn_simulate_PH_path(Lambda, alpha, T_homogeneous, num_phases)  # we generate homogeneous
  trajectory_IPH$times[[i]] <- result$times
  trajectory_IPH$phases[[i]] <- result$phases
  trajectory_IPH$times[[i]] <- sapply(trajectory_IPH$times[[i]], function(x) g_beta_0(x, beta)) # We transform them to non-homogeneous
}
rm(result)

# Discretize paths using delta
for (i in 1:N) {
  result <- fn_discrete_Path(trajectory_IPH$times[[i]], trajectory_IPH$phases[[i]], delta)
  trajectory_IPH$discrete_times[[i]] <- result$times
  trajectory_IPH$discrete_phases[[i]] <- result$phases
}
rm(result)

# Estimation -------------------------------------------------------------------

# Maximum time in the training dataset
max_time_data <- max(sapply(trajectory_IPH$discrete_times, tail, n = 1))  
# Extra time buffer until absorption (10 times the max time in sample)

# Identify chains that reached absorption
trajectory_IPH$absorbed_chain <- sapply(
  trajectory_IPH$discrete_phases, 
  function(phases) tail(phases, 1) == (num_phases + 1)
)

# Get indices for absorbed and non-absorbed chains
absorbed_indices <- which(trajectory_IPH$absorbed_chain)
not_absorbed_indices <- which(!trajectory_IPH$absorbed_chain)

# Initialize Lambda assuming homogeneous times  (i.e., beta_0=0)
EM_results <- fn_estimate_intensity_matrix(
  trajectory_IPH$discrete_times, 
  trajectory_IPH$discrete_phases, 
  num_phases, 
  N
)

Lambda_0 <- EM_results$Lambda_hat
R_vector <- EM_results$R

Lambda_hat <- Lambda_0
Lambda_hat_last <- Lambda_0
Lambda_hat_list[[1]] <- Lambda_0  # Store the initial estimate

# Initialize simulated continuous times and phases (homogeneous assumption)
trajectory_IPH$sim_cont_hom_times <- vector("list", N)
trajectory_IPH$sim_cont_hom_phases <- vector("list", N)

# Process data to generate continuous times and phases (homogeneous assumption)
for (i in seq_len(N)) {
  result <- fn_bridge(Lambda_hat, trajectory_IPH$discrete_phases[[i]], trajectory_IPH$discrete_times[[i]], num_phases) 
  
  # Store simulated times and phases
  trajectory_IPH$sim_cont_hom_times[[i]] <- result$continuous_times
  trajectory_IPH$sim_cont_hom_phases[[i]] <- result$continuous_phases
  
  # Keep ALL available information for non-absorbed chains
  if (i %in% not_absorbed_indices) {
    trajectory_IPH$sim_cont_hom_times[[i]] <- c(trajectory_IPH$sim_cont_hom_times[[i]], tail(trajectory_IPH$discrete_times[[i]], 1))
    trajectory_IPH$sim_cont_hom_phases[[i]] <- c(trajectory_IPH$sim_cont_hom_phases[[i]], tail(trajectory_IPH$discrete_phases[[i]], 1))
  }
}
rm(result)

# Estimate beta_0 using ALL chains (absorbed and non-absorbed). The GD takes care of the non-absorbed ones
estimated_betas <- Gradient_Descent_Complete(
  Lambda = Lambda_hat, 
  beta_0 = beta_0,  
  all_t_times = trajectory_IPH$sim_cont_hom_times, 
  all_phases = trajectory_IPH$sim_cont_hom_phases, 
  eta = eta, 
  error = error_likelihood 
)

# Update beta estimates
beta_0 <- tail(estimated_betas, 1)
beta_hat <- beta_0
beta_hat_list <- beta_hat  # Store estimated betas
number_betas <- length(estimated_betas)

# Transform discrete times to actual "homogeneous times" using the improved beta estimate
trajectory_IPH$transformed_discrete_times <- lapply(
  trajectory_IPH$discrete_times, 
  function(times) sapply(times, g_inv_beta_0, beta_hat)
)


# We improve our estimates for Lambda and beta  
while (iteration <= total_iterations) {
  
  start_time <- Sys.time()  # Measure speed
  
  # Initialize lists to store simulated continuous homogeneous times and phases
  trajectory_IPH$sim_cont_hom_times <- vector("list", N)
  trajectory_IPH$sim_cont_hom_phases <- vector("list", N)
  
  # Bridge by rejection method (using transformed_discrete_times)
  for (i in seq_len(N)){
    result <- fn_bridge(Lambda_hat, trajectory_IPH$discrete_phases[[i]], trajectory_IPH$transformed_discrete_times[[i]], num_phases) 
    trajectory_IPH$sim_cont_hom_times[[i]] <- result$continuous_times
    trajectory_IPH$sim_cont_hom_phases[[i]] <- result$continuous_phases
    
    # Keep ALL available information for non-absorbed chains
    if (i %in% not_absorbed_indices) {
      trajectory_IPH$sim_cont_hom_times[[i]] <- c(trajectory_IPH$sim_cont_hom_times[[i]], tail(trajectory_IPH$transformed_discrete_times[[i]], 1))
      trajectory_IPH$sim_cont_hom_phases[[i]] <- c(trajectory_IPH$sim_cont_hom_phases[[i]], tail(trajectory_IPH$discrete_phases[[i]], 1))
    }
    
  }
  rm(result)
  
  # Simulate absorption times for non-absorbed chains
  if (length(not_absorbed_indices) > 0) {
    for (i in not_absorbed_indices) {
      result <- fn_simulate_absorption_homogenous(Lambda_hat, trajectory_IPH$sim_cont_hom_phases[[i]], trajectory_IPH$sim_cont_hom_times[[i]], num_phases, 10*max_time_data)
      trajectory_IPH$sim_cont_hom_times[[i]] <- result$continuous_times
      trajectory_IPH$sim_cont_hom_phases[[i]] <- result$continuous_phases
    }
    rm(result)
  }

  EM_results <- fn_estimate_intensity_matrix(
    trajectory_IPH$sim_cont_hom_times, 
    trajectory_IPH$sim_cont_hom_phases, 
    num_phases, 
    N
  )
  
  Lambda_hat <- EM_results$Lambda_hat
  R_vector <- EM_results$R
  
  # # Estimate beta
  sim_cont_inhom_times <- lapply(trajectory_IPH$sim_cont_hom_times, function(s_times) {
    g_beta_0(s_times, beta_hat) 
  })
  
  estimated_betas <- Gradient_Descent_Complete(
    Lambda = Lambda_hat, 
    beta_0 = beta_hat, 
    all_t_times = sim_cont_inhom_times, 
    all_phases = trajectory_IPH$sim_cont_hom_phases,
    eta = eta, 
    error = error_likelihood
  )
  
  beta_hat <- tail(estimated_betas, 1)
  
  # Update transformed times with new beta_hat
  trajectory_IPH$transformed_discrete_times <- lapply(trajectory_IPH$discrete_times, 
                                                           function(times) sapply(times, g_inv_beta_0, beta_hat))
  
  elapsed_time <- Sys.time() - start_time
  
  print(Lambda_hat)
  print(elapsed_time)
  cat("Iteration:", iteration, "\n")
  
  Lambda_hat_last <- Lambda_hat
  
  beta_hat_list <- c(beta_hat_list, beta_hat)  # Updates list of estimated betas
  number_betas <- c(number_betas, length(estimated_betas))
  
  Lambda_hat_list[[iteration+1]] <- Lambda_hat  # Updates list of estimated Lambda
  
  iteration <- iteration +1
}
beep(2)

print(beta_hat_list[which(number_betas==2)[1]])
print(beta)
print(Lambda_hat_list[[which(number_betas==2)[1]]])
print(Lambda)

which(number_betas==2)[1]-1


## Plots -------------

plot(beta_hat_list[1:(total_iterations+1)], type = "l", col = "blue", main = "Beta", xlab = "Index", ylab = "Value")
abline(h = beta + 0.01, col = "red", lty = 2, lwd = 2)
abline(h = beta - 0.01, col = "red", lty = 2, lwd = 2)
abline(h = beta, col = "darkgreen", lty = 2, lwd = 2)
abline(h = beta_hat_list[which(number_betas==2)[1]], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[1, 1])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(1,1)", xlab = "Index", ylab = "Value")
abline(h = Lambda[1,1]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[1,1]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[1,1], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][1,1], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[1, 2])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(1,2)", xlab = "Index", ylab = "Value")
abline(h = Lambda[1,2]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[1,2]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[1,2], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][1,2], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[1, 3])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(1,3)", xlab = "Index", ylab = "Value")
abline(h = Lambda[1,3]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[1,3]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[1,3], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][1,3], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[2, 1])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(2,1)", xlab = "Index", ylab = "Value")
abline(h = Lambda[2,1]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[2,1]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[2,1], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][2,1], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[2, 2])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(2,2)", xlab = "Index", ylab = "Value")
abline(h = Lambda[2,2]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[2,2]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[2,2], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][2,2], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[2, 3])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(2,3)", xlab = "Index", ylab = "Value")
abline(h = Lambda[2,3]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[2,3]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[2,3], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][2,3], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[3, 1])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(3,1)", xlab = "Index", ylab = "Value")
abline(h = Lambda[3,1]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[3,1]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[3,1], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][3,1], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[3, 2])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(3,2)", xlab = "Index", ylab = "Value")
abline(h = Lambda[3,2]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[3,2]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[3,2], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][3,2], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[3, 3])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(3,3)", xlab = "Index", ylab = "Value")
abline(h = Lambda[3,3]+0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[3,3]-0.01, col = "red", lty = 2, lwd = 2)
abline(h = Lambda[3,3], col = "darkgreen", lty = 2, lwd = 2)
abline(h =Lambda_hat_list[which(number_betas==2)[1]][[1]][3,3], col = "purple", lty = 2, lwd = 2)


## Improved Plots for Paper

actual_iterations <- length(beta_hat_list[1:(which(number_betas==2)[1])]) # length of the iterations (including initialization)

## --------- beta

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
# ymin <- 0
ymin <- min(c(beta_hat_list[2:actual_iterations], beta))
ymax <- max(c(beta_hat_list[2:actual_iterations], beta))

plot(beta_hat_list[2:actual_iterations],  # We do not include the first element of beta list because it corresponds to initialization
     type = "p", col = "black", pch = 21, bg = "black",
     main = expression(bold("Value of " * hat(beta) * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(beta)),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,  # size of axis labels
     cex.main = 1.5,  # size of main title
     xaxt = "n",
     #yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
#axis(side = 2, at = seq(ymin, ymax, by = 0.2))

abline(h = beta, col = "blue", lty = 1, lwd = 2) # actual value


# ----- Lambda[1,1]

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
lambda_vals <- sapply(Lambda_hat_list, function(x) x[1, 1])[2:(actual_iterations)] # extract the values of Lambda (we do not include the 1st element of list because it corresponds to initialization)
ymin <- min(c(lambda_vals, Lambda[1,1]))
# ymax <- 0
ymax <- max(c(lambda_vals, Lambda[1,1]))


plot(lambda_vals,
     type = "p", col = "black", pch = 21, bg = "black",
     main = expression(bold("Value of " * hat(lambda)[11] * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(lambda)[11]),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,   
     cex.main = 1.5,   
     xaxt = "n",    
     # yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
# axis(side = 2, at = seq(ymin, ymax, by = 0.2))

abline(h = Lambda[1,1], col = "blue", lty = 1, lwd = 2)


# ----- Lambda[1,2]

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
lambda_vals <- sapply(Lambda_hat_list, function(x) x[1, 2])[2:(actual_iterations)] # extract the values of Lambda (we do not include the 1st element of list because it corresponds to initialization)
# ymin <- 0
ymin <- min(c(lambda_vals, Lambda[1,2]))
ymax <- max(c(lambda_vals, Lambda[1,2]))


plot(lambda_vals,
     type = "p", col = "black", pch = 21, bg = "black",
     main = expression(bold("Value of " * hat(lambda)[12] * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(lambda)[12]),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,   
     cex.main = 1.5,   
     xaxt = "n",    
     # yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
# axis(side = 2, at = seq(ymin, ymax, by = 0.1))

abline(h = Lambda[1,2], col = "blue", lty = 1, lwd = 2)


## ------ Lambda[2,1]

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
lambda_vals <- sapply(Lambda_hat_list, function(x) x[2, 1])[2:(actual_iterations)] # extract the values of Lambda (we do not include the 1st element of list because it corresponds to initialization)
# ymin <- 0
ymin <- min(c(lambda_vals, Lambda[2,1]))
ymax <- max(c(lambda_vals, Lambda[2,1]))


plot(lambda_vals,
     type = "p", col = "black", pch = 21, bg = "black",
     main = expression(bold("Value of " * hat(lambda)[21] * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(lambda)[21]),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,   
     cex.main = 1.5,   
     xaxt = "n",    
     # yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
# axis(side = 2, at = seq(ymin, ymax, by = 0.1))

abline(h = Lambda[2,1], col = "blue", lty = 1, lwd = 2)


## ------ Lambda[2,2]

par(mar = c(5, 5, 4, 2))  # expand left margin so the hat (in the label) shows fully
lambda_vals <- sapply(Lambda_hat_list, function(x) x[2, 2])[2:(actual_iterations)] # extract the values of Lambda (we do not include the 1st element of list because it corresponds to initialization)
# ymin <- 0
ymin <- min(c(lambda_vals, Lambda[2,2]))
ymax <- max(c(lambda_vals, Lambda[2,2]))


plot(lambda_vals,
     type = "p", col = "black", pch = 21, bg = "black",
     main = expression(bold("Value of " * hat(lambda)[22] * " through iterations")),
     xlab = expression(paste("Iteration ", i)),
     ylab = expression(hat(lambda)[22]),
     ylim = c(ymin, ymax),
     cex.lab = 1.2,   
     cex.main = 1.5,   
     xaxt = "n",    
     # yaxt = "n"
)

# custom axis ticks
axis(side = 1, at = seq(1, actual_iterations, by = 2))
# axis(side = 2, at = seq(ymin, ymax, by = 0.1))

abline(h = Lambda[2,2], col = "blue", lty = 1, lwd = 2)