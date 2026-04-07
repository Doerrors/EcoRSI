# EcoRSI 
EcoRsi is a package developed to compute the Ecological Restoration Index (RSI) of ecological sites with respect to worst and best sites as reference.


EcoRSI provides a principled, reproducible workflow for computing the **Restoration State Index (RSI)**: a composite, multidimensional measure of how closely a restoration site has converged toward a reference ecosystem, relative to a degraded baseline.

---

## Table of Contents

- [Background](#background)
- [Installation](#installation)
- [How It Works](#how-it-works)
- [Quick Start](#quick-start)
- [Full Worked Example](#full-worked-example)
- [Function Reference](#function-reference)
- [Custom Domain Weights](#custom-domain-weights)
- [Output Structure](#output-structure)
- [Citation](#citation)

---

## Background

Ecological restoration monitoring requires integrating many indicators — vegetation, fauna, soil nutrients, heavy metals — into a single, interpretable measure of recovery progress. EcoRSI addresses this by:

1. Anchoring all indicators to a **Reference Ecosystem (RE)** — what full recovery looks like.
2. Normalizing against a **degraded baseline** (e.g., an Unrestored Mine Dump, UMD) — where recovery started.
3. Grouping indicators into **ecological domains** and computing a weighted composite index.

An **RSI score near 1** means the restoration site closely resembles the reference ecosystem across all measured dimensions. A score **near 0** means it still resembles the degraded baseline.

---

## Installation

EcoRSI is not yet on CRAN. Install directly from GitHub:

```r
# install.packages("devtools")  # if not already installed
devtools::install_github("Doerrors/EcoRSI")
```

Then load it:

```r
library(EcoRSI)
```

---

## How It Works

EcoRSI requires data from at least **three types of sites** measured at the same time points (index):

| Site Type | Role |
|---|---|
| Unrestored Mine Dump (UMD) | Degraded baseline — defines the worst-case bounds |
| Reference Ecosystem (RE) | Target state — all indicators are normalized to this |
| Restoration site(s) | The sites being monitored for recovery |

### Mathematical Framework

**Step 1 — Bounds from the degraded baseline:**

```
L = min(UMD values)   # lower bound, used for positive indicators
U = max(UMD values)   # upper bound, used for negative indicators
```

**Step 2 — Normalize each indicator per time step:**

For *positive* indicators (higher = better, e.g. species richness, soil organic carbon):

```
z⁺ = clamp[0,1]( (X - L) / (RE - L) )
```

For *negative* indicators (lower = better, e.g. heavy metals, salinity):

```
z⁻ = clamp[0,1]( (U - X) / (U - RE) )
```

A score of **1** means the site matches or exceeds the reference. A score of **0** means it is at or below the degraded baseline level.

**Step 3 — Domain sub-indices:**

Indicators are grouped into ecological domains (e.g. Vegetation, Fauna, Soil Nutrients, Heavy Metals). Within each domain, normalized scores are averaged:

```
I_k = mean(z_j)   for all indicators j in domain k
```

**Step 4 — Final RSI:**

```
RSI = Σ (w_k × I_k)   across all K domains
```

By default, equal weights are applied: `w_k = 1/K`. Custom weights can be supplied for context-specific studies.

---

## Quick Start

```r
library(EcoRSI)

# Example: mine dump restoration with 4 ecological domains
result <- computeRSI(
  data                = my_data,
  umd_label           = "UMD",        # degraded baseline label
  re_label            = "RE",         # reference ecosystem label
  domain_list         = list(
    Vegetation  = "floral_richness",
    Fauna       = "faunal_richness",
    SoilNut     = c("N", "P", "K", "OC"),
    HeavyMetal  = c("Mn", "Fe", "pH", "Pb", "Cr")
  ),
  positive_indicators = c("floral_richness", "faunal_richness", "N", "P", "K", "OC"),
  negative_indicators = c("Mn", "Fe", "Pb", "Cr"),
  index_col           = "Index",
  site_col            = "Site"
)

head(result)
```

---

## Full Worked Example

### Step 1 — Prepare your data

Your data must be in **long format**: one row per site per time step, with the same set of time index values repeated across all sites.

```r
# Minimal example dataset
my_data <- data.frame(
  Site            = rep(c("UMD", "RE", "RS1", "RS2"), each = 3),
  Index           = rep(c(1, 5, 10), times = 4),
  floral_richness = c(5, 30, 12, 20,   # UMD, RE, RS1, RS2
                      8, 30, 18, 25,
                      6, 30, 22, 28),
  faunal_richness = c(2, 15, 5,  9,
                      3, 15, 8,  12,
                      2, 15, 11, 14),
  N               = c(0.1, 0.8, 0.3, 0.5,
                      0.1, 0.8, 0.4, 0.6,
                      0.2, 0.8, 0.5, 0.7),
  Pb              = c(120, 10, 80, 50,
                      110, 10, 60, 30,
                      115, 10, 40, 20)
)
```

> **Tip:** Trial datasets are included in the repository (`EcoRSI_trial_dataset_1.xlsx`, etc.). Load them with `readxl::read_excel()`.

---

### Step 2 — Define your domain structure

Group your indicator columns into ecological domains. Every indicator must appear in exactly one domain, and must be declared as either positive or negative.

```r
domain_list <- list(
  Vegetation = "floral_richness",
  Fauna      = "faunal_richness",
  SoilNut    = "N",
  HeavyMetal = "Pb"
)

positive_indicators <- c("floral_richness", "faunal_richness", "N")
negative_indicators <- c("Pb")
```

---

### Step 3 — Run `computeRSI()`

```r
result <- computeRSI(
  data                = my_data,
  umd_label           = "UMD",
  re_label            = "RE",
  domain_list         = domain_list,
  positive_indicators = positive_indicators,
  negative_indicators = negative_indicators,
  index_col           = "Index",
  site_col            = "Site"
)

print(result)
```

---

### Step 4 — Visualize recovery trajectories

```r
library(ggplot2)

# Filter to restoration sites only (exclude UMD and RE)
plot_data <- result[!result$Site %in% c("UMD", "RE"), ]

ggplot(plot_data, aes(x = Index, y = RSI, color = Site, group = Site)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 3) +
  geom_hline(yintercept = 1, linetype = "dashed", color = "darkgreen", linewidth = 0.8) +
  scale_y_continuous(limits = c(0, 1.05)) +
  labs(
    title    = "Ecological Recovery Trajectories (EcoRSI)",
    subtitle = "RSI = 1 indicates full convergence with Reference Ecosystem",
    x        = "Index",
    y        = "Restoration State Index (RSI)",
    color    = "Site"
  ) +
  theme_minimal(base_size = 13)
```

---

## Function Reference

### `computeRSI()`

| Parameter | Type | Description |
|---|---|---|
| `data` | data.frame | Long-format ecological data |
| `umd_label` | character | Label of the degraded baseline site in `site_col` |
| `re_label` | character | Label of the Reference Ecosystem site in `site_col` |
| `domain_list` | named list | Named list mapping domain names to indicator column names |
| `positive_indicators` | character vector | Indicators where higher values = better condition |
| `negative_indicators` | character vector | Indicators where lower values = better condition |
| `index_col` | character | Column name for the time/sampling index |
| `site_col` | character | Column name identifying sites |
| `domain_weights` | named numeric vector | Optional. Custom weights per domain. Must sum to 1. Default: equal weights. |

---

### `validate_EcoRSI_data()`

An internal validation function called automatically by `computeRSI()`. It checks that:

- Both `umd_label` and `re_label` exist in the data
- No missing values are present in the site or index columns
- All sites have an equal number of rows
- The same index values repeat identically across all sites
- No duplicate index values exist within any site

You can also call it directly for pre-flight data checks:

```r
validate_EcoRSI_data(
  data      = my_data,
  site_col  = "Site",
  index_col = "Index",
  umd_label = "UMD",
  re_label  = "RE"
)
# Returns TRUE if all checks pass, stops with a descriptive error otherwise
```

---

## Custom Domain Weights

By default, EcoRSI applies **equal weights** across all domains (`w_k = 1/K`). This reflects the ecological principle that holistic restoration requires simultaneous recovery across all functional domains.

For context-specific studies (e.g. a soil-focused mine restoration), you can supply custom weights:

```r
result_weighted <- computeRSI(
  data                = my_data,
  umd_label           = "UMD",
  re_label            = "RE",
  domain_list         = domain_list,
  positive_indicators = positive_indicators,
  negative_indicators = negative_indicators,
  index_col           = "Index",
  site_col            = "Site",
  domain_weights      = c(
    Vegetation = 0.10,
    Fauna      = 0.10,
    SoilNut    = 0.40,
    HeavyMetal = 0.40
  )
)
```

> ⚠️ **Note:** RSI scores computed with custom weights are **context-specific** and may not be directly comparable with scores from other studies using default equal weights. A warning is issued automatically when custom weights are used.

---

## Output Structure

`computeRSI()` returns a data frame with the following columns:

| Column | Description |
|---|---|
| `Index` | Sampling time index value |
| `Site` | Site label |
| `<DomainName>` | One column per domain — the domain sub-index (I_k), ranging from 0 to 1 |
| `RSI` | Final Restoration State Index — weighted mean of all domain sub-indices |

**Interpretation:**

| RSI | Meaning |
|---|---|
| ~0 | Site condition resembles the degraded baseline |
| 0.5 | Midway between degraded baseline and reference |
| ~1 | Site has converged with the Reference Ecosystem |

---

## Citation

If you use EcoRSI in your research, please cite:

```
Harishwaran V (2025). EcoRSI: Restoration State Index for Ecological Recovery Assessment.
R package. https://github.com/Doerrors/EcoRSI
```

---

## License

This project is open source. See the `DESCRIPTION` file for license details.
