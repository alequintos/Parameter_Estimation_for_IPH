# Gives back a discrete trajectory
# input: Continuous times, continuous phases, Delta,
# Output: Discrete Trajectory: phases and times

fn_discrete_Path=function(cont_times, cont_phases, delta)
{
  times=cont_times
  phases=cont_phases
  n=length(phases)
  t=floor(times[n]/delta)    # Gets the last time and divides it by delta. We round down
  
  discrete_path = matrix(0,(t+2),2)  # We need to add 2 to include 0 and the discrete time after absorption
  discrete_path[,1] = seq(0,t+1)*delta
  discrete_path[1,2] = phases[1]
  
  for(i in 2:(t+1))
  {
    # discrete_path[i,2]=phases[n-length(times[times>(i-1)*delta])]
    discrete_path[i,2]=phases[length(times[times<=(i-1)*delta])]  # It tells you how many of the times are smaller or equal than the current discrete time and then gives you the corresponding phase
  }
  discrete_path[t+2,2]=phases[n]    # The last one should be the exit
  return(list(times=discrete_path[,1], phases=discrete_path[,2]))
}

