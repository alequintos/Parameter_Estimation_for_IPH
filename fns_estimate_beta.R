library(Matrix)

Derivative_Complete_LogLik_Beta <- function(Lambda, beta, all_t_times, all_phases) { # times in the inhomogeneous timeline
  num_phases <- dim(Lambda)[1]
  term1_jumps <- 0      # the derivative of log(h_{\beta}) evaluated at the jumps
  term2_intervals <- 0  # the \lambda_{xx} R_x^{\beta}
  N_paths <- length(all_t_times)
  
  for (j in 1:N_paths) {
    times <- all_t_times[[j]]
    phases <- all_phases[[j]]
    
    # ---------------------------------------------------------
    # TERM 1: Sum over derivative of log(h_{\beta}) evaluated at the jumps
    
    jump_times <- times[-1] 
    
    # Check if the path is absorbed
    is_absorbed <- phases[length(phases)] == (num_phases + 1)
    if(!is_absorbed) {
      jump_times <- jump_times[-length(jump_times)]
    }
    
    if(length(jump_times) > 0) {
      term1_jumps <- term1_jumps + sum(h_beta_0_derivative(jump_times, beta) / h_beta_0(jump_times, beta))
    }
    
    # ---------------------------------------------------------
    # TERM 2: Sum over \lambda_{xx} derivative R_x^{\beta}
    
    for (k in 1:(length(times) - 1)) {
      phase <- phases[k]
      t_start <- times[k]
      t_end <- times[k+1]
      
      # We only accumulate for transient phases
      if (phase <= num_phases) {
        # if t=0, the integral is exactly 0.
        val_start <- ifelse(t_start == 0, 0, int_h_beta_0_derivative(t_start, beta))
        val_end   <- int_h_beta_0_derivative(t_end, beta)
        
        d_R_beta <- val_end - val_start
        term2_intervals <- term2_intervals + Lambda[phase, phase] * d_R_beta
      }
    }
  }
  
  return(term1_jumps + term2_intervals)
}


# ---------------------------------------

Gradient_Descent_Complete <- function(Lambda, beta_0, all_t_times, all_phases, eta, error) {
  num_phases <- dim(Lambda)[1]
  difference_likelihoods <- Inf
  betas <- beta_0
  
  m <- length(betas)
  
  # Calculate initial proportional log-likelihood
  current_ll <- Complete_LogLik_Beta_Proportional(Lambda, beta_0, all_t_times, all_phases)
  
  while(difference_likelihoods > error) {
    last_beta <- tail(betas, 1)
    print(last_beta)
    
    # 1. Compute the gradient 
    grad <- Derivative_Complete_LogLik_Beta(Lambda, last_beta, all_t_times, all_phases)
    
    # 2. Update beta
    new_beta <- max(10^(-5), last_beta + eta * grad)
    
    # 3. Compute the new proportional log-likelihood
    new_ll <- Complete_LogLik_Beta_Proportional(Lambda, new_beta, all_t_times, all_phases)
    
    # 4. Update the difference and lists
    difference_likelihoods <- abs(new_ll - current_ll)
    current_ll <- new_ll
    betas <- c(betas, new_beta)
    
    m <- length(betas)
    
    print(tail(betas,1))
  }
  
  return(betas)
}

# --------------------------


Complete_LogLik_Beta_Proportional <- function(Lambda, beta, all_t_times, all_phases) {
  num_phases <- dim(Lambda)[1]
  term1_jumps <- 0
  term2_intervals <- 0
  N_paths <- length(all_t_times)
  
  for (j in 1:N_paths) {
    times <- all_t_times[[j]]
    phases <- all_phases[[j]]
    
    # TERM 1: Sum over (log h(t)) evaluated at the jumps
    jump_times <- times[-1] 
    
    # Check if the path is absorbed
    is_absorbed <- phases[length(phases)] == (num_phases + 1)
    if(!is_absorbed) {
      jump_times <- jump_times[-length(jump_times)]
    }
    
    if(length(jump_times) > 0) {
      term1_jumps <- term1_jumps + sum(log(h_beta_0(jump_times, beta)))
    }
    
    # TERM 2: Sum over Intervals (Integral of h(t)) 
    
    for (k in 1:(length(times) - 1)) {
      phase <- phases[k]
      t_start <- times[k]
      t_end <- times[k+1]
      
      # We only accumulate for transient phases
      if (phase <= num_phases) {
        
        # Safe integration: if t=0, the integral is exactly 0.
        val_start <- ifelse(t_start == 0, 0, int_h_beta_0(t_start, beta))
        val_end   <- int_h_beta_0(t_end, beta)
        
        d_R_beta <- val_end - val_start
        term2_intervals <- term2_intervals + Lambda[phase, phase] * d_R_beta
      }
    }
  }
  
  return(term1_jumps + term2_intervals)
}