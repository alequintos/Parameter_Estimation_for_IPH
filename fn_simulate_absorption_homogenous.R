# Simulates what is missing of a continuous path until it reaches absorption


# Q = Approximation for subintensity matrix, phases = phases at continuous times
# times = times at continuous times, num_phases = number of phases of subintensity matrix
# simulation_time = simulation time. This is not too important because we will repeat until the chain is absorbed

fn_simulate_absorption_homogenous <- function(Q, phases, times, num_phases, simulation_time) {
  # getting information to start simulation
  terminal_phase <- phases[length(phases)]
  terminal_time <- times[length(times)]
  
  # we create a vector to start of the simulation at the last phase
  alpha_simulation <- numeric(num_phases)
  alpha_simulation[terminal_phase] <- 1
  
  is_absorbed <- 0  # we keep track if chain is absorbed or not
  iteration <- 1
  
  # We simulate until absorption
  while(is_absorbed == 0){
    result <- fn_simulate_PH_path(Q, alpha_simulation, iteration*simulation_time, num_phases)
    
    if ( tail(result$phases,1) == (num_phases+1) ){  # The chain was absorbed
      is_absorbed <- 1
    }
    iteration <- iteration + 1
  }
  
  # We put together what we had with what we simulated
  continuous_times <- c(times, terminal_time + result$times[-1])   # we do not add the first position of our simulation
  continuous_phases <- c(phases, result$phases[-1])
  
  list(continuous_times = continuous_times, continuous_phases = continuous_phases)
}