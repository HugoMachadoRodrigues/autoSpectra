# Compute the Q-statistic for each row of a PCA residual matrix

Q_i = e_i %\*% e_i^T (eq. 2 in Safanelli 2025).

## Usage

``` r
q_statistic(E)
```

## Arguments

- E:

  Numeric matrix of PCA residuals (X - X_hat)

## Value

Numeric vector of Q values, length nrow(E)
