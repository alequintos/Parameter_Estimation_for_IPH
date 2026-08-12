# Auxiliary functions

# Yslas and real data ----------------------

h_beta_0 <- function(t, beta) { exp(beta*t) }
h_beta_0_derivative <- function(t, beta) { t * exp(beta*t) }

int_h_beta_0 <- function(t, beta) {
  (exp(beta * t) - 1) / beta
}

int_h_beta_0_derivative <- function(t, beta) {
  (exp(beta * t)*(beta*t-1) + 1) / beta^2
}

# Cumulative intensity function
g_inv_beta_0 <- function(t, beta) {
  (exp(beta*t)-1)/beta
}

# Inverse cumulative intensity function
g_beta_0 <- function(g_inv_value, beta) {  # (Yslas paper)
  log(beta*g_inv_value + 1)/beta
}


# Bladt and Peralta paper ----------------------

# h_beta_0 <- function(t, beta) { beta * t^(beta-1)}
# h_beta_0_derivative <- function(t, beta) {t^(beta-1)*(1 + beta*log(t))}
# 
# # Integrals of h and derivative of h wrt beta
# int_h_beta_0 <- function(t, beta) {
#   t^beta
# }
# 
# int_h_beta_0_derivative <- function(t, beta) {
#   t^beta * log(t)
# }
# 
# # Cumulative intensity function
# g_inv_beta_0 <- function(t, beta) {
#   t^beta
# }
# 
# # Inverse cumulative intensity function
# g_beta_0 <- function(g_inv_value, beta) {
#   g_inv_value^(1/beta)
# }
