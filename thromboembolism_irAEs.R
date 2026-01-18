


library(dplyr)
library(readr)
library(tidyr)
library(stringr)

# ----------------------------
# 0) Read input files
# ----------------------------
pt_counts    <- read_csv("RSTUDIOdiseases.csv", show_col_types = FALSE)
drugs_counts <- read_csv("RSTUDIOdrugs.csv",    show_col_types = FALSE)

# Make column names safe (e.g., spaces -> dots)
names(pt_counts)    <- make.names(names(pt_counts))
names(drugs_counts) <- make.names(names(drugs_counts))


# ----------------------------
# 1)  nij, ni, nj, n
# ----------------------------
# In a standard 2x2 table:
#   nij = a
#   ni  = a + b
#   nj  = a + c
#   n   = a + b + c + d


# ----------------------------
# 2) Metrics: PRR + PRR CI, IC + IC025
# ----------------------------
# PRR uses the full 2x2 table (a,b,c,d).
# IC (Information Component, BCPNN-style) uses nij and expected count E = ni*nj/n.
# We include a 0.5 continuity correction (common in practice) to stabilize small counts.
calc_metrics <- function(a, b, c, d, nij, ni, nj, n) {
  # --- PRR ---
  prr <- (a / (a + b)) / (c / (c + d))
  se_log_prr <- sqrt(1/a - 1/(a + b) + 1/c - 1/(c + d))
  prr_l <- exp(log(prr) - 1.96 * se_log_prr)
  prr_u <- exp(log(prr) + 1.96 * se_log_prr)
  
  # --- IC ---
  # Expected count under independence:
  #   E = (ni * nj) / n
  E <- (ni * nj) / n
  
  # Continuity correction:
  nij_cc <- nij + 0.5
  E_cc   <- E + 0.5
  
  # Information Component:
  #   IC = log2( (nij_cc) / (E_cc) )
  IC <- log2(nij_cc / E_cc)
  
  # Approximate standard error for IC (normal approximation on log scale):
  # var(log(nij/E)) ~ 1/nij + 1/E
  # convert natural log variance to log2 variance by dividing by (ln 2)^2
  var_IC <- (1 / nij_cc + 1 / E_cc) / (log(2)^2)
  se_IC  <- sqrt(var_IC)
  
  IC025 <- IC - 1.96 * se_IC
  
  tibble(
    PRR = prr, PRR025 = prr_l, PRR975 = prr_u,
    IC = IC, IC025 = IC025,
    Expected = E
  )
}

# ----------------------------
# 3) Table 2 PT-level
# ----------------------------
table2_rep <- pt_counts %>%
  mutate(
    a   = as.numeric(a),
    b   = as.numeric(b),
    c   = as.numeric(c),
    d   = as.numeric(d),
    nij = as.numeric(nij),
    ni  = as.numeric(ni),
    nj  = as.numeric(nj),
    n   = as.numeric(n)
  ) %>%
  rowwise() %>%
  mutate(
    metrics = list(
      calc_metrics(a, b, c, d, nij, ni, nj, n)
    )
  ) %>%
  unnest(metrics) %>%
  ungroup() %>%
  select(
    PT,
    a, b, c, d,
    nij, ni, nj, n,
    Expected,
    PRR, PRR025, PRR975,
    IC, IC025
  )

write_csv(table2_rep, "Table2_reproduction.csv")

# ----------------------------
# 4) Table 3 Drug x PT
# ----------------------------
table3_rep <- drugs_counts %>%
  mutate(
    a   = as.numeric(a),
    b   = as.numeric(b),
    c   = as.numeric(c),
    d   = as.numeric(d),
    nij = as.numeric(nij),
    ni  = as.numeric(ni),
    nj  = as.numeric(nj),
    n   = as.numeric(n)
  ) %>%
  rowwise() %>%
  mutate(
    metrics = list(
      calc_metrics(a, b, c, d, nij, ni, nj, n)
    )
  ) %>%
  unnest(metrics) %>%
  ungroup() %>%
  select(
    drugs,
    PT,
    a, b, c, d,
    nij, ni, nj, n,
    Expected,
    PRR, PRR025, PRR975,
    IC, IC025
  )

write_csv(table3_rep, "Table3_reproduction.csv")
