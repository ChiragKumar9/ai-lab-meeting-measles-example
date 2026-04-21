# =============================================================================
# Age-Stratified Measles Model with Waning Vaccine Immunity
# =============================================================================
#
# Purpose:
#   Template script for assessing whether measles vaccine immunity wanes over
#   time using an age-stratified compartmental (ODE) model.
#
# Model structure: MSEIR + V (vaccinated with waning immunity)
#
# Compartments per age group:
#   M  - Maternally immune (youngest age group only)
#   S  - Susceptible
#   E  - Exposed / latent (infected but not yet infectious)
#   I  - Infectious
#   R  - Recovered (permanently immune via natural infection)
#   V  - Vaccinated (immune, but subject to waning)
#
# Key transitions:
#   Births -> M  (fraction with maternal immunity)
#   Births -> S  (fraction without maternal immunity)
#   M -> S       : loss of maternal immunity (rate delta)
#   S -> E       : infection driven by age-structured contact matrix
#   E -> I       : progression through latent period (rate sigma)
#   I -> R       : recovery (rate gamma)
#   S -> V       : vaccination (age-specific rate * vaccine efficacy)
#   V -> S       : waning vaccine immunity (rate omega) -- key parameter
#   All compartments subject to age-specific mortality and aging
#
# Usage:
#   1. Adjust parameters in Section 2 to match your setting.
#   2. Optionally replace the placeholder contact matrix with an empirical
#      matrix (e.g. from Mossong et al. 2008, POLYMOD study).
#   3. Run the script to simulate and visualise dynamics.
#   4. Section 8 compares multiple waning-rate scenarios.
#
# Dependencies: deSolve, ggplot2, dplyr, tidyr
#   Install with: install.packages(c("deSolve","ggplot2","dplyr","tidyr"))
#
# References:
#   Anderson & May (1991) Infectious Diseases of Humans. Oxford University Press.
#   McLean & Blower (1993) Imperfect vaccines and the emergence of measles.
#     Proc R Soc B 253:9-13.
#   Mossong et al. (2008) Social contacts and mixing patterns relevant to the
#     spread of infectious diseases. PLoS Med 5(3):e74.
# =============================================================================

library(deSolve)
library(ggplot2)
library(dplyr)
library(tidyr)

# =============================================================================
# 1. Age Group Configuration
# =============================================================================

age_group_labels <- c(
  "<1 yr", "1-4 yr", "5-9 yr", "10-14 yr",
  "15-19 yr", "20-29 yr", "30-49 yr", "50+ yr"
)
n_age_groups <- length(age_group_labels)

# =============================================================================
# 2. Model Parameters
# =============================================================================
# Time unit: days throughout (convert annual rates where noted).

params <- list(

  n_age_groups = n_age_groups,

  # ---- Demography --------------------------------------------------------

  # Total population size
  n_total = 1e6,

  # Fraction of population in each age group (must sum to 1)
  age_fractions = c(0.012, 0.050, 0.065, 0.065, 0.065, 0.130, 0.260, 0.353),

  # Per-capita birth rate (per day)
  birth_rate = 14 / 1000 / 365, # ~14 per 1000 per year

  # Age-specific daily mortality rates (per capita per day)
  death_rates = c(
    5.0e-5, # <1 yr
    2.0e-6, # 1-4 yr
    5.0e-7, # 5-9 yr
    5.0e-7, # 10-14 yr
    8.0e-7, # 15-19 yr
    1.5e-6, # 20-29 yr
    5.0e-6, # 30-49 yr
    3.0e-5  # 50+ yr
  ),

  # Rate of aging out of each group (per day) = 1 / (width of age class in days)
  # The last group has no aging out (rate = 0).
  aging_rates = c(
    1 / (1 * 365),   # <1 yr   -> 1-4 yr  (class width: 1 yr)
    1 / (4 * 365),   # 1-4 yr  -> 5-9 yr  (class width: 4 yr)
    1 / (5 * 365),   # 5-9 yr  -> 10-14 yr
    1 / (5 * 365),   # 10-14 yr -> 15-19 yr
    1 / (5 * 365),   # 15-19 yr -> 20-29 yr
    1 / (10 * 365),  # 20-29 yr -> 30-49 yr
    1 / (20 * 365),  # 30-49 yr -> 50+ yr
    0                # 50+ yr  (open-ended; no aging out)
  ),

  # ---- Disease Natural History -------------------------------------------

  # Rate of progression from E to I: 1 / latent period
  # Measles latent period ~12 days
  sigma = 1 / 12,

  # Rate of recovery from I to R: 1 / infectious period
  # Measles infectious period ~8 days
  gamma = 1 / 8,

  # ---- Transmission ------------------------------------------------------

  # Baseline transmission rate (per infectious contact per day).
  # Scaled by the contact matrix.  Tune so that R0 ~15 in an unvaccinated pop.
  beta = 0.5,

  # Contact matrix (n_age_groups x n_age_groups).
  # Entry [i, j] = mean daily contacts a person in age group i has with
  # people in age group j.
  # NULL here -- populated below by build_contact_matrix().
  # Replace with an empirical matrix (e.g. POLYMOD) for a real analysis.
  contact_matrix = NULL,

  # ---- Vaccination -------------------------------------------------------

  # Age-specific vaccination rates (per susceptible per day).
  # These represent the rate at which susceptibles in each age class
  # receive an effective vaccine dose.
  vaccination_rates = c(
    0,               # <1 yr      (too young for MMR)
    1 / (1 * 365),   # 1-4 yr     (primary MMR at ~12 months -> ~1-yr class)
    1 / (10 * 365),  # 5-9 yr     (second dose / catch-up)
    1 / (20 * 365),  # 10-14 yr
    1 / (50 * 365),  # 15-19 yr
    1 / (100 * 365), # 20-29 yr
    1 / (200 * 365), # 30-49 yr
    1 / (500 * 365)  # 50+ yr
  ),

  # Probability that a vaccine dose successfully confers immunity
  vaccine_efficacy = 0.97,

  # ---- Waning Vaccine Immunity (key parameter) ---------------------------

  # Rate of loss of vaccine-induced immunity (per day).
  # omega = 1 / (duration of vaccine protection in days).
  # Set omega = 0 for lifelong (non-waning) immunity.
  # Example values:
  #   1 / (30*365) -> protection wanes over ~30 years
  #   1 / (15*365) -> protection wanes over ~15 years
  #   1 / (10*365) -> protection wanes over ~10 years
  omega = 1 / (30 * 365), # default: slow waning over 30 years

  # ---- Maternal Immunity -------------------------------------------------

  # Rate of loss of maternally-derived immunity (per day)
  # Duration of maternal immunity ~6 months
  delta = 1 / 180,

  # Fraction of newborns that receive maternal immunity
  # (depends on population-level immunity of mothers)
  prop_immune_mothers = 0.90
)

# -----------------------------------------------------------------------------
# Build a simplified assortative contact matrix
# (replace with an empirical matrix for a real analysis)
# -----------------------------------------------------------------------------

build_contact_matrix <- function(n) {
  # Baseline off-diagonal contacts (between age groups)
  cm <- matrix(0.5, nrow = n, ncol = n)
  # Higher within-group contacts (assortative mixing)
  diag(cm) <- c(2, 8, 6, 5, 4, 3, 3, 2)[seq_len(n)]
  # Symmetrise (reciprocal contacts)
  cm <- (cm + t(cm)) / 2
  cm
}

params$contact_matrix <- build_contact_matrix(n_age_groups)

# =============================================================================
# 3. ODE Model Function
# =============================================================================
#
# State vector layout (length = 6 * n_age_groups):
#   M[1:n], S[n+1:2n], E[2n+1:3n], I[3n+1:4n], R[4n+1:5n], V[5n+1:6n]

measles_ode <- function(t, state, params) {
  with(params, {
    n <- n_age_groups

    # Unpack state
    M <- state[seq_len(n)]
    S <- state[(n + 1):(2 * n)]
    E <- state[(2 * n + 1):(3 * n)]
    I <- state[(3 * n + 1):(4 * n)]
    R <- state[(4 * n + 1):(5 * n)]
    V <- state[(5 * n + 1):(6 * n)]

    N <- M + S + E + I + R + V # total population per age group

    # Force of infection: lambda[i] = beta * sum_j C[i,j] * I[j] / N[j]
    # Avoid division by zero in empty age groups
    prevalence <- I / ifelse(N > 0, N, 1)
    lambda <- beta * as.vector(contact_matrix %*% prevalence)

    # Initialise derivative vectors
    dM <- numeric(n)
    dS <- numeric(n)
    dE <- numeric(n)
    dI <- numeric(n)
    dR <- numeric(n)
    dV <- numeric(n)

    total_population <- sum(N)
    births_per_day <- birth_rate * total_population

    for (i in seq_len(n)) {

      # --- Aging flows ---
      # Inflow from younger age group (i-1); none for the youngest group (i=1)
      aging_in_M <- if (i > 1) aging_rates[i - 1] * M[i - 1] else 0
      aging_in_S <- if (i > 1) aging_rates[i - 1] * S[i - 1] else 0
      aging_in_E <- if (i > 1) aging_rates[i - 1] * E[i - 1] else 0
      aging_in_I <- if (i > 1) aging_rates[i - 1] * I[i - 1] else 0
      aging_in_R <- if (i > 1) aging_rates[i - 1] * R[i - 1] else 0
      aging_in_V <- if (i > 1) aging_rates[i - 1] * V[i - 1] else 0

      # --- Birth inflows (youngest group only) ---
      birth_M <- if (i == 1) prop_immune_mothers * births_per_day else 0
      birth_S <- if (i == 1) (1 - prop_immune_mothers) * births_per_day else 0

      # --- Vaccination flow (S -> V) ---
      vaccinations <- vaccination_rates[i] * vaccine_efficacy * S[i]

      # --- ODEs ---

      # Maternally immune
      dM[i] <- birth_M +
        aging_in_M -
        delta * M[i] -         # loss of maternal immunity -> S
        death_rates[i] * M[i] -
        aging_rates[i] * M[i]

      # Susceptible
      dS[i] <- birth_S +
        aging_in_S +
        delta * M[i] +         # from M
        omega * V[i] -         # waning vaccine immunity -> S (key waning term)
        lambda[i] * S[i] -     # infection -> E
        vaccinations -         # vaccination -> V
        death_rates[i] * S[i] -
        aging_rates[i] * S[i]

      # Exposed (latent)
      dE[i] <- aging_in_E +
        lambda[i] * S[i] -     # new infections from S
        sigma * E[i] -         # progression -> I
        death_rates[i] * E[i] -
        aging_rates[i] * E[i]

      # Infectious
      dI[i] <- aging_in_I +
        sigma * E[i] -         # progression from E
        gamma * I[i] -         # recovery -> R
        death_rates[i] * I[i] -
        aging_rates[i] * I[i]

      # Recovered (permanently immune via natural infection)
      dR[i] <- aging_in_R +
        gamma * I[i] -         # recovery from I
        death_rates[i] * R[i] -
        aging_rates[i] * R[i]

      # Vaccinated (immune, waning)
      dV[i] <- aging_in_V +
        vaccinations -         # newly vaccinated from S
        omega * V[i] -         # waning immunity -> S
        death_rates[i] * V[i] -
        aging_rates[i] * V[i]
    }

    list(c(dM, dS, dE, dI, dR, dV))
  })
}

# =============================================================================
# 4. Initial Conditions
# =============================================================================

set_initial_conditions <- function(params) {
  n <- params$n_age_groups
  n0 <- params$n_total * params$age_fractions

  # Approximate starting immune fractions by age
  # (older cohorts more likely to be naturally immune or previously vaccinated)
  prop_natural_immune <- c(0.00, 0.05, 0.20, 0.40, 0.60, 0.70, 0.80, 0.85)
  prop_vaccinated <- c(0.00, 0.70, 0.85, 0.88, 0.85, 0.80, 0.75, 0.70)

  # Seed infectious individuals in the 5-9 yr age group
  seed_I <- numeric(n)
  seed_I[3] <- 5 # 5 infectious in 5-9 yr group
  seed_E <- seed_I * 2

  r0 <- n0 * prop_natural_immune
  v0 <- n0 * prop_vaccinated * (1 - prop_natural_immune)
  m0 <- c(n0[1] * params$prop_immune_mothers, numeric(n - 1))
  i0 <- seed_I
  e0 <- seed_E
  s0 <- pmax(n0 - m0 - r0 - v0 - i0 - e0, 0)

  state0 <- c(m0, s0, e0, i0, r0, v0)

  comp_names <- c("M", "S", "E", "I", "R", "V")
  names(state0) <- unlist(lapply(
    comp_names,
    function(comp) paste0(comp, "_", params$age_group_labels)
  ))

  state0
}

params$age_group_labels <- age_group_labels
state0 <- set_initial_conditions(params)

# =============================================================================
# 5. Run Baseline Simulation
# =============================================================================

# Simulate 30 years (in days)
t_end <- 30 * 365
times <- seq(0, t_end, by = 1)

cat("Running baseline simulation...\n")
sol_baseline <- ode(
  y = state0,
  times = times,
  func = measles_ode,
  parms = params,
  method = "lsoda"
)

sol_df <- as.data.frame(sol_baseline)

# =============================================================================
# 6. Process Output
# =============================================================================

# Helper: reshape to long format for one compartment
extract_compartment <- function(df, comp, age_labels) {
  cols <- paste0(comp, "_", age_labels)
  df_long <- df |>
    select(all_of(c("time", cols))) |>
    pivot_longer(
      cols = all_of(cols),
      names_to = "age_group",
      values_to = "count"
    ) |>
    mutate(
      age_group = factor(
        sub(paste0(comp, "_"), "", age_group),
        levels = age_labels
      ),
      compartment = comp
    )
  df_long
}

comp_names <- c("M", "S", "E", "I", "R", "V")
results <- bind_rows(lapply(
  comp_names,
  extract_compartment,
  df = sol_df,
  age_labels = age_group_labels
))

# Total infectious over time
total_I_df <- data.frame(
  time_years = sol_df$time / 365,
  total_I = rowSums(sol_df[, grep("^I_", names(sol_df))])
)

# =============================================================================
# 7. Visualisation
# =============================================================================

# -- 7a. Total infectious prevalence over time --
p_total <- ggplot(total_I_df, aes(x = time_years, y = total_I)) +
  geom_line(color = "firebrick", linewidth = 0.8) +
  labs(
    title = "Measles: Total Infectious Individuals Over Time",
    subtitle = paste0("Waning rate omega = ", round(params$omega * 365, 4), " /yr"),
    x = "Time (years)",
    y = "Number Infectious"
  ) +
  theme_bw()

print(p_total)

# -- 7b. Infectious by age group --
p_age_I <- results |>
  filter(compartment == "I") |>
  ggplot(aes(x = time / 365, y = count, color = age_group)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Measles: Infectious Individuals by Age Group",
    x = "Time (years)",
    y = "Number Infectious",
    color = "Age group"
  ) +
  theme_bw()

print(p_age_I)

# -- 7c. Vaccinated fraction by age group (shows waning) --
p_vax <- results |>
  filter(compartment %in% c("V", "S", "E", "I", "R", "M")) |>
  group_by(time, age_group) |>
  summarise(
    V = sum(count[compartment == "V"]),
    N = sum(count),
    .groups = "drop"
  ) |>
  mutate(vax_fraction = V / ifelse(N > 0, N, 1)) |>
  ggplot(aes(x = time / 365, y = vax_fraction, color = age_group)) +
  geom_line(linewidth = 0.7) +
  labs(
    title = "Vaccinated Fraction by Age Group (Waning Immunity)",
    x = "Time (years)",
    y = "Fraction Vaccinated",
    color = "Age group"
  ) +
  theme_bw()

print(p_vax)

# =============================================================================
# 8. Waning Immunity Scenarios
# =============================================================================
# Compare total measles burden under different vaccine waning assumptions.

waning_scenarios <- list(
  "No waning (lifelong)"    = 0,
  "Slow waning (~30 yr)"    = 1 / (30 * 365),
  "Moderate waning (~15 yr)" = 1 / (15 * 365),
  "Fast waning (~10 yr)"    = 1 / (10 * 365)
)

cat("Running waning immunity scenario comparisons...\n")

run_scenario <- function(omega_val, scenario_name, params, state0, times) {
  p <- params
  p$omega <- omega_val
  sol <- ode(
    y = state0,
    times = times,
    func = measles_ode,
    parms = p,
    method = "lsoda"
  )
  sol_df <- as.data.frame(sol)
  data.frame(
    time_years = sol_df$time / 365,
    total_I = rowSums(sol_df[, grep("^I_", names(sol_df))]),
    scenario = scenario_name,
    stringsAsFactors = FALSE
  )
}

scenario_df <- bind_rows(
  mapply(
    run_scenario,
    omega_val = waning_scenarios,
    scenario_name = names(waning_scenarios),
    MoreArgs = list(params = params, state0 = state0, times = times),
    SIMPLIFY = FALSE
  )
)

p_scenarios <- ggplot(
  scenario_df,
  aes(x = time_years, y = total_I, color = scenario)
) +
  geom_line(linewidth = 0.8) +
  labs(
    title = "Effect of Waning Vaccine Immunity on Measles Dynamics",
    x = "Time (years)",
    y = "Total Infectious",
    color = "Waning scenario"
  ) +
  theme_bw() +
  theme(legend.position = "bottom")

print(p_scenarios)

# =============================================================================
# 9. Summary Statistics
# =============================================================================

cat("\n--- Baseline simulation summary ---\n")
cat("Total population (initial):", round(sum(state0)), "\n")
cat("Initial infectious:", round(sum(state0[grep("^I_", names(state0))])), "\n")
cat("Peak total infectious:", round(max(total_I_df$total_I)), "\n")
cat("Time to peak (years):", round(total_I_df$time_years[which.max(total_I_df$total_I)], 2), "\n")

scenario_peaks <- scenario_df |>
  group_by(scenario) |>
  summarise(
    peak_I = max(total_I),
    time_to_peak_yr = time_years[which.max(total_I)],
    .groups = "drop"
  )

cat("\n--- Peak infectious by waning scenario ---\n")
print(scenario_peaks)
