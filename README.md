# Age-Stratified Measles Model with Waning Vaccine Immunity

Template R script for assessing whether measles vaccine immunity wanes over time
using an age-stratified compartmental (ODE) model.

## Model Overview

The model is an **MSEIR + V** system with **8 age groups**
(`<1 yr`, `1-4 yr`, `5-9 yr`, `10-14 yr`, `15-19 yr`, `20-29 yr`, `30-49 yr`, `50+ yr`).

Compartments per age group:

| Compartment | Description |
|-------------|-------------|
| **M** | Maternally immune (youngest group only) |
| **S** | Susceptible |
| **E** | Exposed / latent |
| **I** | Infectious |
| **R** | Recovered (permanently immune via natural infection) |
| **V** | Vaccinated (immune, subject to waning) |

The key parameter for this analysis is **`omega`** — the rate at which
vaccine-induced immunity wanes (`V → S`).  Section 8 of the script
compares four waning scenarios ranging from lifelong immunity (`omega = 0`)
to immunity that fades over ~10 years.

Age-structured transmission is driven by a contact matrix (a simplified
assortative-mixing matrix is provided; replace with an empirical matrix such
as [POLYMOD](https://doi.org/10.1371/journal.pmed.0050074) for a real analysis).

## How to Run

1. **Install dependencies** (one-time setup):

   ```r
   install.packages(c("deSolve", "ggplot2", "dplyr", "tidyr"))
   ```

2. **Run the model script**:

   ```r
   source("model_measles_waning.R")
   ```

   Or open `model_measles_waning.R` in RStudio and run it interactively.

3. **Customise**:
   - Adjust epidemiological parameters in **Section 2** (e.g. `omega`, `beta`,
     `vaccination_rates`).
   - Replace the placeholder contact matrix with an empirical one
     (e.g. from Mossong et al. 2008).
   - Add or modify waning scenarios in **Section 8**.

## Outputs

- **Plot 1** – Total infectious prevalence over time (baseline scenario).
- **Plot 2** – Infectious individuals by age group.
- **Plot 3** – Vaccinated fraction by age group, illustrating waning over time.
- **Plot 4** – Comparison of measles dynamics under different waning-rate
  assumptions.
- **Console summary** – Peak infectious count and time-to-peak for each
  waning scenario.

## References

- Anderson & May (1991) *Infectious Diseases of Humans*. Oxford University Press.
- McLean & Blower (1993) Imperfect vaccines and the emergence of measles.
  *Proc R Soc B* 253:9-13.
- Mossong et al. (2008) Social contacts and mixing patterns relevant to the
  spread of infectious diseases. *PLoS Med* 5(3):e74.
