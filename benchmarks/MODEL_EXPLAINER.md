# DLNM Model Explainer — Beginner's Guide

A step-by-step walkthrough of the dataset, every model input, and how the model
gets constructed, grounded in the actual code and data.

---

## The Dataset

The original dataset is `chicagoNMMAPS` — daily records for the city of Chicago
from 1987 to 2000. **One row = one day.**

```
        date  time  year  month  doy       dow   death   temp    pm10
  1987-01-01     1  1987      1    1  Thursday    130   -0.28   26.96
  1987-01-02     2  1987      1    2    Friday    150    0.56      NA
  1987-01-03     3  1987      1    3  Saturday    101    0.56   32.84
  ...
  2000-12-31  5114  2000     12  366    Sunday     97   -2.78   18.21
```

The columns we care about for the model:

| Column  | What it is |
|---------|-----------|
| `death` | Total deaths in Chicago on that day. This is what we're trying to explain — the **outcome**. |
| `temp`  | Average temperature in °C that day — the **exposure** we're studying. |
| `time`  | Just a counter: day 1, 2, 3 … 5114. Used to track the long-run trend. |
| `dow`   | Day of the week (Monday, Tuesday, …). People die at different rates on weekends vs weekdays. |

The benchmark datasets are just this table replicated across 24, 235, or 2350
fake "cities" with slight random differences in temperature and death rates, with
a `city` column added.

---

## The Scientific Question

**Does temperature affect how many people die — and if so, does that effect play
out immediately or over several days?**

For example: if it's extremely hot today, do people die today? Or mostly 3 days
later? Or spread over 2 weeks? The DLNM is designed to answer exactly this.

---

## Building the Model — Step by Step

The model formula in the benchmark is:

```r
glm(death ~ cb + ns(time, df = time_df) + dow,
    family = quasipoisson())
```

There are three ingredients on the right-hand side.

---

### Ingredient 1: `cb` — the Cross-Basis (the DLNM part)

This is the core of the package. It gets built in two sub-steps before the GLM
ever runs.

#### Step 1a: The Variable Basis — "what shape is the temperature effect?"

Temperature isn't measured in a vacuum. A 1°C increase from 5°C to 6°C might
have a different effect than a 1°C increase from 32°C to 33°C. The **variable
basis** captures this non-linear shape by transforming the raw temperature column
into a set of smooth curve columns.

For config C2, we use `ns(df=5)` — a **natural spline with 5 degrees of
freedom**. Think of it like fitting the data with 5 overlapping curved tent
functions that together can approximate any smooth curve shape. The single `temp`
column becomes 5 new columns, each representing how "activated" that part of the
temperature curve is on a given day.

```
  raw temp →  basisvar (5 columns)
  -0.28    →  [0.12, 0.45, 0.31, 0.09, 0.02]
   0.56    →  [0.11, 0.44, 0.32, 0.10, 0.02]
  33.33    →  [0.00, 0.01, 0.08, 0.42, 0.61]
```

#### Step 1b: The Lag Basis — "how does the effect decay over time?"

The effect of temperature on death doesn't only happen on the same day. If it's
hot today, some people might die today, some tomorrow, some next week. The **lag
basis** models this decay shape.

For C2, we use a lag range of 0–21 days and `ns(df=4)` — another natural spline,
but this time applied to the lag dimension (0, 1, 2, … 21). It transforms the 22
possible lag days into 4 smooth columns that can represent any decay pattern:
sharp and immediate, gradual, or U-shaped.

#### Step 1c: The Tensor Product — combining both

`crossbasis()` combines the two bases via a **tensor product**: for every
combination of (variable basis column, lag basis column), it computes one output
column. C2: 5 × 4 = **20 columns**.

Each of those 20 columns answers the question: "given the lagged temperature
history, how much of exposure pattern V with lag pattern L is present today?"

Concretely, for each day `i`, the code looks back 21 days, takes the temperature
value each of those days, weights it by the lag basis, and produces the
cross-basis row for day `i`. That's what the nested loop in `crossbasis.R` does:

```r
# For each of the 5 variable basis columns:
#   Build a 22-column lag matrix by shifting the basis column by 0,1,2...21 days
#   Then multiply by each of the 4 lag basis columns
#   → produces 5 × 4 = 20 output columns
for v in 1..5:
    lag_matrix = Lag(basisvar[, v], 0:21)   # shape: (n_days × 22)
    for l in 1..4:
        cb[, 4*(v-1)+l] = lag_matrix %*% basislag[, l]
```

This is the **primary bottleneck** identified in benchmarking. At 1GB (12M rows),
`Lag()` takes ~49 seconds just to build those shifted matrices, accounting for
~75% of total `crossbasis()` time. It is the top target for Rust/SIMD
optimisation (P1 in the benchmark plan).

The end result `cb` is a matrix with one row per day and 20 columns. Each row
encodes the entire 21-day temperature history of that day in a compressed, smooth
form.

---

### Ingredient 2: `ns(time, df = time_df)` — the Long-Run Trend

`time` is just the day counter (1 through 5114). Over 14 years, death rates
slowly drift up and down due to things that have nothing to do with temperature:
aging population, flu seasons, medical advances, etc.

If we don't account for this, the model might mistake a long-run seasonal pattern
for a temperature effect.

`ns(time, df = time_df)` fits a natural spline over the trend in time. `time_df`
controls how flexible that trend curve is allowed to be:

| `time_df` | Flexibility | Model matrix memory at 1GB |
|-----------|------------|---------------------------|
| 7         | One bend per year — broad seasonal/annual trends only | ~1.7 GB |
| 14        | One bend per 6 months | ~2.4 GB |
| 98 (= 7 × 14 years) | One bend per ~7 weeks — standard epidemiological practice | ~10.5 GB (OOM at 1GB scale) |

**This is what `time_df` is.** It's how many "knots" (bends) the long-run time
trend curve is allowed to have. Epidemiologists typically use 7 per year (one per
~7 weeks) to fully capture seasonality — that's where 98 comes from (7 × 14
years of data). For the 1GB benchmarks, it was reduced to 7 total to keep the
model matrix within R's 32 GB vector memory limit.

---

### Ingredient 3: `dow` — Day of Week

People behave differently on weekends: hospitals have fewer staff, people are
less likely to seek care, etc. A simple factor variable with 7 levels
(Sunday–Saturday). R automatically converts it into 6 binary columns (one level
is the reference category).

---

### Putting it all Together: the GLM

The full model matrix for one day looks like this (C2, `time_df = 7`):

```
[ cb_v1l1, cb_v1l2, ..., cb_v5l4  |  ns(time)_1...7  |  Mon, Tue, Wed, Thu, Fri, Sat ]
[      20 cross-basis columns      |    7 trend cols   |      6 day-of-week cols       ]
= 33 columns total per row
```

`glm(..., family = quasipoisson())` then fits this by **Iteratively Reweighted
Least Squares (IRLS)**:

1. Start with a guess for all 33 coefficients
2. Compute predicted death counts from the current coefficients
3. Compute how far off the predictions are (residuals)
4. Re-weight each row by how uncertain the prediction is
5. Solve a weighted least-squares problem to update the coefficients
6. Repeat until convergence (4 iterations in the 1GB benchmarks)

The result is 33 fitted coefficients — 20 of which belong to the cross-basis and
encode the full 2D temperature–lag–mortality surface.

---

### After the GLM: `crosspred()` and `crossreduce()`

Once we have the 20 cross-basis coefficients, `crosspred()` reconstructs the
interpretable 2D surface:

- "At 33°C vs 15°C, what is the relative risk of death at lag 0? At lag 3? At
  lag 10?"
- "Accumulated over the full 21-day lag window, what is the overall excess risk
  of a hot day?"

`crossreduce()` then collapses that 2D surface down to 1D — either the overall
cumulative effect across all lags, or the effect at one specific lag across all
temperatures.

These two steps are fast (< 1 second at any scale) because they operate only on
the 20 fitted coefficients and their uncertainty estimates — not on the raw data.

---

## Model Configurations Used in Benchmarking

| Config | Variable basis | Lag basis | Lag range | CB columns | Purpose |
|--------|---------------|-----------|-----------|------------|---------|
| C1 | `lin` (1 col) | `poly(degree=4)` (5 cols) | 0–15 days | 5 | Minimal / linear DLM |
| C2 | `ns(df=5)` (5 cols) | `ns(df=4)` (4 cols) | 0–21 days | 20 | Typical epidemiology |
| C3 | `bs(df=6)` (6 cols) | `ns(df=4)` (4 cols) | 0–40 days | 24 | Extended lag window |
| C4 | `ps(df=10)` (10 cols) | `ps(df=5)` (5 cols) | 0–30 days | 50 | Penalised high-df |
| C5 | `ps(df=15)` (15 cols) | `ps(df=8)` (8 cols) | 0–60 days | 120 | Stress test |

More cross-basis columns = richer model, but also more memory and compute at
every stage.

---

## Benchmark Results Summary

### End-to-end stage times (median seconds)

| Scale | Rows | Config | `crossbasis` | `glm` | `crosspred` | `crossreduce` |
|-------|------|--------|-------------|-------|------------|--------------|
| 10 MB | 122K | C1 | 0.03s | 3.66s | 0.02s | 0.02s |
| 10 MB | 122K | C2 | 0.22s | 4.43s | 0.06s | 0.05s |
| 100 MB | 1.2M | C1 | 0.38s | 14.6s | 0.10s | 0.11s |
| 100 MB | 1.2M | C2 | 2.15s | 17.7s | 0.10s | 0.11s |
| 1 GB | 12M | C1 | 5.49s | 35.3s† | ~0.5s* | ~0.8s* |
| 1 GB | 12M | C2 | 64.6s | 74.2s† | ~0.2s* | ~0.2s* |

\* Extrapolated via power-law scaling from 10MB and 100MB results.
† GLM run with `time_df=7` instead of the standard 98 — standard formula
  exceeds R's 32 GB vector memory limit at this scale (model matrix alone is
  ~10.5 GB for C1, ~11.9 GB for C2).

### `crossbasis` internal hotspots at 1GB

| Config | `lag_matrix` | `multiply_loop` | Total |
|--------|-------------|----------------|-------|
| C1 (5 cols) | 4.08s (74%) | 0.60s (11%) | 5.49s |
| C2 (20 cols) | 48.6s (75%) | 7.94s (12%) | 64.6s |

Lag matrix construction dominates at every scale and is the primary target for
the planned Rust/SIMD rewrite.
