# Load required libraries
library(odin)
library(deSolve)
library(ggplot2)
library(dplyr)

# Define the SIR model using odin's domain-specific language
sir_model <- odin({
  # Derivatives
  deriv(S) <- -beta * S * I / N
  deriv(I) <- beta * S * I / N - gamma * I
  deriv(R) <- gamma * I
  
  # Initial conditions
  initial(S) <- S_ini
  initial(I) <- I_ini
  initial(R) <- R_ini
  
  # Parameters
  beta <- user(0.5)     # transmission rate
  gamma <- user(0.1)    # recovery rate
  
  # Total population (constant)
  N <- S_ini + I_ini + R_ini
  
  # Initial values (user-defined)
  S_ini <- user(999)
  I_ini <- user(1)
  R_ini <- user(0)
  
  # Output variables
  output(prevalence) <- I / N
  output(attack_rate) <- R / N
})

# Create model instance with parameters
model <- sir_model(
  beta = 0.3,      # transmission rate
  gamma = 0.1,     # recovery rate (1/infectious period)
  S_ini = 999,     # initial susceptible
  I_ini = 1,       # initial infected
  R_ini = 0        # initial recovered
)

# Set time points for simulation
times <- seq(0, 100, by = 0.1)

# Run the simulation
result <- model$run(times)

# Convert to data frame for plotting
sir_data <- as.data.frame(result)

# Create visualization
plot_sir <- ggplot(sir_data, aes(x = t)) +
  geom_line(aes(y = S, color = "Susceptible"), size = 1) +
  geom_line(aes(y = I, color = "Infected"), size = 1) +
  geom_line(aes(y = R, color = "Recovered"), size = 1) +
  scale_color_manual(values = c("Susceptible" = "blue", 
                               "Infected" = "red", 
                               "Recovered" = "green")) +
  labs(title = "SIR Model Dynamics",
       x = "Time",
       y = "Number of Individuals",
       color = "Compartment") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(plot_sir)

# Plot prevalence over time
plot_prevalence <- ggplot(sir_data, aes(x = t, y = prevalence)) +
  geom_line(color = "red", size = 1) +
  labs(title = "Disease Prevalence Over Time",
       x = "Time",
       y = "Prevalence (I/N)") +
  theme_minimal()

print(plot_prevalence)

# Calculate basic reproduction number (R0)
R0 <- 0.3 / 0.1  # beta / gamma
cat("Basic Reproduction Number (R0):", R0, "\n")

# Summary statistics
peak_infected <- max(sir_data$I)
peak_time <- sir_data$t[which.max(sir_data$I)]
final_attack_rate <- tail(sir_data$attack_rate, 1)

cat("Peak number of infected:", peak_infected, "\n")
cat("Time of peak infection:", peak_time, "\n")
cat("Final attack rate:", round(final_attack_rate * 100, 2), "%\n")

# Function to run model with different parameters
run_sir_scenario <- function(beta, gamma, S_ini = 999, I_ini = 1, R_ini = 0, max_time = 100) {
  model <- sir_model(beta = beta, gamma = gamma, 
                     S_ini = S_ini, I_ini = I_ini, R_ini = R_ini)
  times <- seq(0, max_time, by = 0.1)
  result <- model$run(times)
  return(as.data.frame(result))
}

# Example: Compare different transmission rates
scenarios <- data.frame(
  beta = c(0.2, 0.3, 0.4),
  gamma = c(0.1, 0.1, 0.1),
  scenario = c("Low transmission", "Medium transmission", "High transmission")
)

# Run scenarios
scenario_results <- list()
for(i in 1:nrow(scenarios)) {
  result <- run_sir_scenario(scenarios$beta[i], scenarios$gamma[i])
  result$scenario <- scenarios$scenario[i]
  scenario_results[[i]] <- result
}

# Combine results
all_scenarios <- do.call(rbind, scenario_results)

# Plot comparison
comparison_plot <- ggplot(all_scenarios, aes(x = t, y = I, color = scenario)) +
  geom_line(size = 1) +
  labs(title = "SIR Model: Effect of Different Transmission Rates",
       x = "Time",
       y = "Number of Infected",
       color = "Scenario") +
  theme_minimal() +
  theme(legend.position = "bottom")

print(comparison_plot)