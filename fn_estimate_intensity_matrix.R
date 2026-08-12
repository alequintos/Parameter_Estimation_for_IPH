fn_estimate_intensity_matrix <- function(all_times, all_phases, num_phases, N){
  # We estimate the intensity matrix
  
  # Initialize a list to store the time spent in each phase for each path
  time_spent_in_phases <- vector("list", N)
  
  for (i in 1:N) {
    # Calculate the time differences between phase changes
    time_differences <- diff(all_times[[i]])
    
    # Store the time spent in each phase
    time_spent_in_phases[[i]] <- data.frame(
      phase = all_phases[[i]][-length(all_phases[[i]])], # Exclude the last phase as it is the absorbing state
      time_spent = time_differences
    )
  }
  
  time_in_phase <- vector("list", num_phases)
  R <- numeric(num_phases)
  
  for (i in 1:num_phases){
    # Initialize a list to store the time spent in phase i across all paths (We need a list because not every path will visit the state i or some paths will visit more than once)
    time_in_phase[[i]]<- list()
    for (j in 1:N){
      # Filter the data frame to get only the rows where the phase is i
      phase_k_times <- time_spent_in_phases[[j]]$time_spent[time_spent_in_phases[[j]]$phase == i]
      
      # Append these times to the list
      time_in_phase[[i]] <- c(time_in_phase[[i]], phase_k_times)
    }
    # Convert the list to a numeric vector
    time_in_phase[[i]] <- unlist(time_in_phase[[i]])
    
    R[i] <- sum(time_in_phase[[i]])
  }
  
  
  
  # Initialize a matrix to store transition counts between phases
  N_matrix <- matrix(0, nrow = num_phases, ncol = num_phases+1)
  
  # Loop through each path
  for (i in 1:N) {
    # Get the sequence of phases for the i-th path
    phases <- all_phases[[i]]
    
    # Loop through consecutive phase pairs
    for (j in 1:(length(phases) - 1)) {
      phase_k <- phases[j]
      phase_l <- phases[j + 1]
      
      # Increment the count for the transition from phase k to phase l
      N_matrix[phase_k, phase_l] <- N_matrix[phase_k, phase_l] + 1
    }
  }
  
  estimated_Lambda <- matrix(0, nrow = num_phases, ncol = num_phases+1)
  for (i in 1:num_phases){
    for (j in 1:(num_phases+1)){
      if (i != j){
        estimated_Lambda[i,j] <- N_matrix[i,j]/R[i]
      }
    }
    estimated_Lambda[i,i]<- - sum(estimated_Lambda[i, ])
  }
  
  # Return a list containing the matrix, the occupation times (R), and the transition counts (N_matrix) 
  return(list(
    Lambda_hat = estimated_Lambda[1:num_phases, 1:num_phases],
    R = R,
    N_matrix = N_matrix
  ))
  
}