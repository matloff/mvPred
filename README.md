
# mvPred: A Central Repo for Handling Missing Values in Predictive Applications

Methods for handling missing values abound, but they are almost all
focused in applications in which the goal is effect estimation, e.g.
calculation and statistical inference for regression coefficients.
Our focus here is instead on prediction.

Note that this should not be confused with packages that use regression
methods for inputation. There, the missing values in a variable
X<sub>i</sub> may be inputed using predictions in which X<sub>2</sub>
is regressed on the other X<sub>j</sub>. Instead, we are interested in
applications in prediction itself is the focus, such as forecasting or
disease diagnosis.

## Installation

You can install `mvPred` from CRAN with:

```r
install.packages("mvPred")
```

Then load the package with:

```r
library(mvPred)
```

## Main Methods

The package currently includes methods for predictive modeling with missing data, including:

- complete-case analysis via `bootstrap(..., method = "CC")`
- available-case regression via `bootstrap(..., method = "AC")`
- PREFILL-based approaches via `bootstrap(..., method = "PREFILL")`
- TOWER-based modeling via `bootstrap(..., method = "TOWER")`

For PREFILL, the supported imputation methods include:

- `mice`
- `Amelia`
- `missForest`
- complete-case imputation through `impute_method = "complete"`

Core modeling helpers in the package include:

- `bootstrap()`
- `lm_ac()`
- `lm_prefill()`
- `lm_tower()`

## Example

```r
library(mvPred)

data("auto-mpg", package = "mvPred")

df <- auto_mpg
df$car_name <- NULL

for (nm in names(df)) {
  suppressWarnings(df[[nm]] <- as.numeric(df[[nm]]))
}

res <- bootstrap(
  data = df,
  yName = "mpg",
  k = 5,
  task = "regression",
  method = "CC"
)

res$RMSE_mean
```

## Package Data

The package also includes datasets used in examples and experiments, such as:

- `auto-mpg`
- `english`
- `NHkids`


