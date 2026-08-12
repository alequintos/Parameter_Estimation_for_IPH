# Simulates a continuous path when we have a starting and end point

# This function allows for various lengths of bridges by incorporating a dynamic delta. 
# It is particularly useful when we are working with inhomogeneous cases

# Q = Approximation for subintensity matrix, discrete_phases = phases at discrete times
# discrete_times = times at discrete times, num_phases = number of phases of subintensity matrix

fn_bridge <- function(Q, discrete_phases, discrete_times, num_phases) {
  number_discrete_times <- length(discrete_phases)    # number of times we will repeat the procedure
  
  # They will store the continuous version of the MJP. We initialize them with 0 (starting time) and the first given phase
  continuous_times <- c(0)
  continuous_phases <- discrete_phases[1]
  
  dynamic_delta <- diff(discrete_times)  # Stores the length of the bridges
  
  j <- 1
  while (j < number_discrete_times){
    initial_phase <- discrete_phases[j]
    terminal_phase <- discrete_phases[j+1]
    
    vector_initial_phase <- numeric(num_phases)  # Creates an "artificial vector" to start the MJP in every iteration
    vector_initial_phase[initial_phase] <- 1
    
    result <- fn_simulate_PH_path(Q, vector_initial_phase, dynamic_delta[j], num_phases) # We simulate the process in the given time interval (dynamic_delta[j])
    
    end<-length(result$phases)
    
    if (result$phases[end]==terminal_phase){      # We will accept the simulation 
      if (end >1){                                # There was an actual jump and so, we store its value
        continuous_times <- c(continuous_times, discrete_times[j]+result$times[-1])   # we do not add the first position of our simulation
        continuous_phases <- c(continuous_phases, result$phases[-1])
      }
      j <- j+1
    }
    else{
    }
  }
  list(continuous_times = continuous_times, continuous_phases = continuous_phases)
}