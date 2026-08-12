# Generates sample paths necessary for the estimation.
# It generates the paths at continuous times and then discretizes according to delta

source("fn_simulate_PH_path.R")
source("fns_auxiliary.R")     # For g and g_inv

source("fn_discrete_path.R")   # to discretize

fn_generate_simulated_data_IPH <- function(alpha, Lambda, beta, N, T_inhomogeneous, delta) {
  T_homogeneous <- g_inv_beta_0(T_inhomogeneous, beta)  # Transforms time to homogenous
  
  # Number of transient phases (excluding the absorbing state)
  num_phases <- length(alpha)
  
  # Initialize storage for trajectories
  trajectory_IPH <- vector("list", 4)
  trajectory_IPH[[1]] <- vector("list", N)
  trajectory_IPH[[2]] <- vector("list", N)
  trajectory_IPH[[3]] <- vector("list", N)
  trajectory_IPH[[4]] <- vector("list", N)
  names(trajectory_IPH) <- c("times", "phases", "discrete_times", "discrete_phases")
  
  # Generate sample paths and transform them to inhomogeneous
  for (i in 1:N) {
    result <- fn_simulate_PH_path(Lambda, alpha, T_homogeneous, num_phases)
    trajectory_IPH$times[[i]] <- sapply(result$times, function(x) g_beta_0(x, beta)) # Transform times
    trajectory_IPH$phases[[i]] <- result$phases
  }
  
  # Discretize the paths
  for (i in 1:N) {
    result <- fn_discrete_Path(trajectory_IPH$times[[i]], trajectory_IPH$phases[[i]], delta)
    trajectory_IPH$discrete_times[[i]] <- result$times
    trajectory_IPH$discrete_phases[[i]] <- result$phases
  }
  
  # Return the full trajectory object. It's a list
  return(trajectory_IPH)
}