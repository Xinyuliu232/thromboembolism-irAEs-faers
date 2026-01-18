# Reproduction of FAERS Disproportionality Analysis for Thromboembolism irAEs

This repository documents the reproduction of statistical analyses conducted in my **previously completed pharmacovigilance study** on thromboembolism-related immune-related adverse events (irAEs) associated with immune checkpoint inhibitors (ICIs).

The original analysis was performed several years ago.  
Due to time lapse and incomplete preservation of intermediate files, this project was created to **reconstruct, document, and validate the analytical workflow** in a fully reproducible manner.

---

## Background

The original study investigated thromboembolic adverse events associated with immune checkpoint inhibitors using the FDA Adverse Event Reporting System (FAERS) database.

The primary statistical methods included:

- Proportional Reporting Ratio (PRR)
- Information Component (IC, Bayesian Confidence Propagation Neural Network)

This repository does **not aim to replicate another group’s work**, but rather to **reproduce my own prior analysis** to ensure methodological transparency, reproducibility, and long-term documentation.

---

## Data description

Two pre-aggregated contingency tables are used as input:

- `RSTUDIOdiseases.csv`  
  PT-level contingency tables for overall ICI exposure

- `RSTUDIOdrugs.csv`  
  Drug × PT contingency tables for individual immune checkpoint inhibitors

Each table includes the following counts:

- `a`, `b`, `c`, `d` (2×2 contingency table)
- `nij`, `ni`, `nj`, `n` (observed and marginal counts)

These tables represent the finalized analytical inputs used in the original study.

---

## Note on early-stage data processing

During the original project period, parts of the early-stage data summarization and tabulation were performed for exploratory and descriptive purposes.

In this reproduction, all statistical metrics are recalculated programmatically in R to ensure:

- transparency
- reproducibility
- long-term maintainability

The focus of this repository is therefore **methodological reconstruction rather than raw database reprocessing**.

---

## Statistical methods

For each drug–event combination, disproportionality metrics were calculated as follows:

### Expected count

\[
E_{ij} = \frac{n_i \times n_j}{n}
\]

### Proportional Reporting Ratio (PRR)

\[
PRR = \frac{a/(a+b)}{c/(c+d)}
\]

with 95% confidence intervals calculated on the log scale.

### Information Component (IC)

\[
IC = \log_2\left(\frac{n_{ij} + 0.5}{E_{ij} + 0.5}\right)
\]

A continuity correction of 0.5 was applied to stabilize small counts.

The lower bound of the 95% credibility interval (IC025) was estimated using a normal approximation.

---

## How to run

```r
source("thromboembolism_irAEs.R")
