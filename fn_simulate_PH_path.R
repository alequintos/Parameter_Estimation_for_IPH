# Function to simulate one sample path (it returns times the process jumps and the phases visited)

# Lambda = subintensity matrix (excluding absorbing state), alpha = initial probabilities
# T = maximum simulation time, num_phases = Number of phases (excluding absorption)

fn_simulate_PH_path <- function(Lambda, alpha, T, num_phases) {
  current_phase <- which(rmultinom(1, 1, alpha) == 1) # Start in phase based on initial distribution
  current_time <- 0
  times <- c(0)
  phases <- current_phase
  
  next_phase <- 0
  
  repeat {
    # Calculate the holding time in the current phase
    lambda <- -Lambda[current_phase, current_phase] # Rate of leaving current phase
    if(lambda == 0) break # No transitions out of this phase
    
    holding_time <- rexp(1, lambda) # Exponential holding time
    current_time <- current_time + holding_time # Update cumulative time
    
    if(current_time >= T) break
    
    # Determine the next phase, including the possibility of absorption
    transition_probabilities <- c(Lambda[current_phase,], -sum(Lambda[current_phase,])) # Last term (the added one) is absorption
    transition_probabilities[current_phase] <- 0 # No transition to the same phase
    transition_probabilities <- transition_probabilities / sum(transition_probabilities) # Normalize probabilities
    
    # To protect from an error when estimating
    if(transition_probabilities[num_phases+1]<0 && abs(transition_probabilities[num_phases+1]) < 1e-10){  # Only if the probability is negative and small
      transition_probabilities[num_phases+1]<- 0
    }
    
    
    next_phase <- which(rmultinom(1, 1, transition_probabilities) == 1) # Determine next phase
    
    if(next_phase == num_phases + 1) {
      # Transition to the absorbing state
      break
    } else {
      current_phase <- next_phase # Transition to the next transient phase
      times <- c(times, current_time)
      phases <- c(phases, current_phase)
    }
  }
  
  # Include the final absorbing state if reached
  if (next_phase == num_phases + 1) {
    times <- c(times, current_time)
    phases <- c(phases, num_phases + 1) # Absorbing state
  }
  
  list(times = times, phases = phases)  # Function returns times the process jumps and the phases visited
}
