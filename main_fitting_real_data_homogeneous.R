# Set seed for reproducibility
set.seed(2)

# Importing and cleaning data --------------------

# install.packages("msm")

# Load necessary library and load data
library(msm)
data_raw <- cav  

# Extract unique identifiers
identifiers <- data_raw[[1]]

# Create a named list
data_example <- setNames(
  lapply(data_raw[-1], function(element) split(element, identifiers)),
  names(data_raw)[-1]
)

# Remove unnecessary objects
rm(data_raw, identifiers)


# Estimation ------------------

# Load required libraries ------------------------------------------------
# install.packages("beepr")  # To make notification sounds

library(expm)   # For matrix exponentials
library(beepr)  # To make notification sounds

# Source necessary scripts ------------------------------------------------
source("fn_estimate_intensity_matrix.R")         # Estimate lambda
source("fn_bridge.R")     
source("fn_simulate_absorption_homogenous.R")
source("fn_simulate_PH_path.R")


# Estimation parameters ---------------------------------------------------
alpha <- c(1, 0, 0)  # Initial state probabilities

total_iterations <- 300  # number of iterations in estimation algorithm 

num_phases <- length(alpha)
elapsed_time <- 0
iteration <- 1
Lambda_hat_list <- list()


# Process data ------------------------------------------------------------
# Identify absorbed chains in original data
absorbed_chain <- sapply(data_example$state, function(x) tail(x, 1) == (num_phases+1))

# Indices of absorbed and non-absorbed chains
absorbed_indices_original_data <- which(absorbed_chain)
not_absorbed_indices_original_data <- which(!absorbed_chain)  # considers all indices

# Filter non-absorbed indices to include only those visiting all three states (1,2,3)
contains_all_states <- sapply(data_example$state[not_absorbed_indices_original_data], function(x) all(c(1:num_phases) %in% x))
valid_not_absorbed_indices <- not_absorbed_indices_original_data[contains_all_states]

# Number of absorbed and non-absorbed chains
N_absorbed <- length(absorbed_indices_original_data)
N_not_absorbed <- length(valid_not_absorbed_indices)

# Select half of absorbed and non-absorbed indices for training
train_absorbed <- sample(absorbed_indices_original_data, size = floor(N_absorbed/2), replace = FALSE)
train_non_absorbed <- sample(valid_not_absorbed_indices, size = floor(N_not_absorbed/2), replace = FALSE)
indices_training_data <- c(train_absorbed, train_non_absorbed)

# Remaining indices for testing
indices_testing_data <- setdiff(c(absorbed_indices_original_data, valid_not_absorbed_indices), indices_training_data)

# Remove patients with duplicated absorption times: the two-sample KS test
# requires continuous data, so tied observations are dropped before testing.
# See Section 5 of the paper.
indices_training_data <- indices_training_data[!indices_training_data %in% c(499, 525)]

# Remove patients with duplicated absorption times: the two-sample KS test
# requires continuous data, so tied observations are dropped before testing.
# See Section 5 of the paper.
indices_testing_data <- indices_testing_data[!indices_testing_data %in% c(129, 215, 584, 621)]


# Dataset sizes
size_data <- length(data_example$years)
size_training_data <- length(indices_training_data)
size_testing_data <- length(indices_testing_data)


# Prepare training data for estimation ------------------------------------
data_for_estimation <- list(
  discrete_times = data_example$years[indices_training_data],
  discrete_phases = data_example$state[indices_training_data]
)

# Estimation -------------------------------------------------------------------

# Maximum time in the training dataset
max_time_data <- max(sapply(data_for_estimation$discrete_times, tail, n = 1))  
# Extra time buffer until absorption (10 times the max time in sample)

# Identify chains that reached absorption
data_for_estimation$absorbed_chain <- sapply(
  data_for_estimation$discrete_phases, 
  function(phases) tail(phases, 1) == (num_phases + 1)
)

# Get indices for absorbed and non-absorbed chains
absorbed_indices <- which(data_for_estimation$absorbed_chain)
not_absorbed_indices <- which(!data_for_estimation$absorbed_chain)

# Initialize Lambda
Lambda_0 <- fn_estimate_intensity_matrix(data_for_estimation$discrete_times, data_for_estimation$discrete_phases, num_phases, size_training_data)$Lambda_hat
Lambda_hat <- Lambda_0
Lambda_hat_last <- Lambda_0
Lambda_hat_list[[1]] <- Lambda_0  # Store the initial estimate

# We improve our estimates for Lambda  
while (iteration <= total_iterations) {
  
  start_time <- Sys.time()  # Measure speed
  
  # Initialize lists to store simulated continuous homogeneous times and phases
  data_for_estimation$sim_cont_hom_times <- vector("list", size_training_data)
  data_for_estimation$sim_cont_hom_phases <- vector("list", size_training_data)
  
  # Bridge by rejection method
  for (i in seq_len(size_training_data)){
    result <- fn_bridge(Lambda_hat, data_for_estimation$discrete_phases[[i]], data_for_estimation$discrete_times[[i]], num_phases) 
    data_for_estimation$sim_cont_hom_times[[i]] <- result$continuous_times
    data_for_estimation$sim_cont_hom_phases[[i]] <- result$continuous_phases
    
    # Adjust times for absorbed chains (the real data has actual absorption times)
    if (i %in% absorbed_indices) {
      data_for_estimation$sim_cont_hom_times[[i]][length(data_for_estimation$sim_cont_hom_times[[i]])] <- 
        tail(data_for_estimation$discrete_times[[i]], 1)
    }
    
    # Keep ALL available information for non-absorbed chains
    if (i %in% not_absorbed_indices) {
      data_for_estimation$sim_cont_hom_times[[i]] <- c(data_for_estimation$sim_cont_hom_times[[i]], tail(data_for_estimation$discrete_times[[i]], 1))
      data_for_estimation$sim_cont_hom_phases[[i]] <- c(data_for_estimation$sim_cont_hom_phases[[i]], tail(data_for_estimation$discrete_phases[[i]], 1))
    }
    
  }
  rm(result)
  
  # Simulate absorption times for non-absorbed chains
  if (length(not_absorbed_indices) > 0) {
    for (i in not_absorbed_indices) {
      result <- fn_simulate_absorption_homogenous(Lambda_hat, data_for_estimation$sim_cont_hom_phases[[i]], data_for_estimation$sim_cont_hom_times[[i]], num_phases, 10*max_time_data)
      data_for_estimation$sim_cont_hom_times[[i]] <- result$continuous_times
      data_for_estimation$sim_cont_hom_phases[[i]] <- result$continuous_phases
    }
    rm(result)
  }
  
  # Estimate Lambda using the homogeneous times
  Lambda_hat <- fn_estimate_intensity_matrix(data_for_estimation$sim_cont_hom_times, 
                                             data_for_estimation$sim_cont_hom_phases, 
                                             num_phases, 
                                             size_training_data)$Lambda_hat
  
  elapsed_time <- Sys.time() - start_time
  
  print(Lambda_hat)
  print(elapsed_time)
  cat("Iteration:", iteration, "\n")
  
  Lambda_hat_last <- Lambda_hat
  
  Lambda_hat_list[[iteration+1]] <- Lambda_hat  # Updates list of estimated Lambda
  
  iteration <- iteration +1
}
beep(2)



# Our proposal for point estimators
print(Reduce("+", Lambda_hat_list[(total_iterations-19):(total_iterations+1)])/length(Lambda_hat_list[(total_iterations-19):(total_iterations+1)]))

# Comparing ----------------------------------------------------------------
alpha2 <- alpha
Lambda2 <- Reduce("+", Lambda_hat_list[(total_iterations-19):(total_iterations+1)])/length(Lambda_hat_list[(total_iterations-19):(total_iterations+1)])

num_phases <- length(alpha)
num_simulations <- length(intersect(absorbed_indices_original_data, indices_testing_data))

# Initialize list for storing simulations
trajectory_IPH2 <- list(times = vector("list", num_simulations),
                        phases = vector("list", num_simulations))

for(i in 1:num_simulations) {
  result <- fn_simulate_PH_path(Lambda2, alpha2, 100000, num_phases)  # we generate homogeneous
  trajectory_IPH2$times[[i]] <- result$times
  trajectory_IPH2$phases[[i]] <- result$phases
}
rm(result)

sample_IPH2 <- sapply(trajectory_IPH2$times, function(x) tail(x, 1))

sample_IPH_test <-  sapply(intersect(absorbed_indices_original_data, indices_testing_data), function(i) tail(data_example$years[[i]], 1))
# sample_IPH_train <-  sapply(intersect(absorbed_indices_original_data, indices_training_data), function(i) tail(data_example$years[[i]], 1))
sample_IPH_test <- unique(sample_IPH_test)

# Histogram Plot -------------------------------------------------
hist(sample_IPH_test, probability = TRUE, col = rgb(0, 0, 1, 0.5), border = "black", breaks = "fd",
     main = "Comparison of Real vs Estimated Distribution", xlab = "Time", ylab = "Density")

hist(sample_IPH2, probability = TRUE, col = rgb(1, 0, 0, 0.5), border = "black", breaks = "fd", add = TRUE)

legend("topright", legend = c("Real Data", "Simulated Data"), fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), bty = "n")


hist(sample_IPH2, probability = TRUE, col = rgb(1, 0, 0, 0.5), border = "black", breaks = "fd",
     main = "Comparison of Real vs Estimated Distribution", xlab = "Time", ylab = "Density")
hist(sample_IPH_test, probability = TRUE, col = rgb(0, 0, 1, 0.5), border = "black", breaks = "fd", add = TRUE)

legend("topright", legend = c("Real Data", "Simulated Data"), fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), bty = "n")

# Kolmogorov-Smirnov Test -------------------------------------------------
ks_test_result <- ks.test(sample_IPH_test, sample_IPH2)
print(ks_test_result)

if (ks_test_result$p.value > 0.05) {
  message("KS test: We cannot reject the null hypothesis -> distributions are similar.")
} else {
  message("KS test: Significant difference detected (p < 0.05).")
}



plot(sapply(Lambda_hat_list, function(x) x[1, 1])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(1,1)", xlab = "Index", ylab = "Value")
abline(h = Lambda2[1,1], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[1, 2])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(1,2)", xlab = "Index", ylab = "Value")
abline(h = Lambda2[1,2], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[2, 1])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(2,1)", xlab = "Index", ylab = "Value")
abline(h = Lambda2[2,1], col = "purple", lty = 2, lwd = 2)

plot(sapply(Lambda_hat_list, function(x) x[2, 2])[1:(total_iterations+1)], type = "l", col = "blue", main = "Lambda(2,2)", xlab = "Index", ylab = "Value")
abline(h = Lambda2[2,2], col = "purple", lty = 2, lwd = 2)


## Improved Histograms (for paper) ------

# pdf("Rplot_RealData_Homogeneous.pdf",
#     width = 570/72,
#     height = 326/72,
#     pointsize = 12)

# Compute combined x-axis limits
x_limits <- range(c(sample_IPH_test, sample_IPH2))

# Compute histograms without plotting (to get max density)
h1 <- hist(sample_IPH_test, breaks = "fd", plot = FALSE)
h2 <- hist(sample_IPH2, breaks = "fd", plot = FALSE)

# Get the maximum density from both
y_max <- max(c(h1$density, h2$density))

# Now plot with ylim that accommodates both
hist(sample_IPH_test, probability = TRUE, col = rgb(0, 0, 1, 0.5), border = "black", breaks = "fd",
     main = "Real vs. Simulated Absorption Time Distributions",
     xlab = "Time", ylab = "Density", xlim = x_limits, ylim = c(0, y_max))

# Add second histogram
hist(sample_IPH2, probability = TRUE, col = rgb(1, 0, 0, 0.5), border = "black", breaks = "fd", add = TRUE)

# Add legend
legend("topright", legend = c("Real Data", "Simulated Data"),
       fill = c(rgb(0, 0, 1, 0.5), rgb(1, 0, 0, 0.5)), bty = "n")

# dev.off()

