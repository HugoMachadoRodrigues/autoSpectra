# Upper Confidence Limit for the Q-statistic (Jackson-Mudholkar 1979)

Implements eq. 3 of Safanelli 2025: \$\$UCL_Q = \theta_1 \left\[1 -
\theta_2 h_0 \frac{1 - h_0}{\theta_1^2} + \frac{z\_\alpha (2 \theta_2
h_0^2)^{1/2}}{\theta_1}\right\]^{1/h_0}\$\$ where \\\theta_k =
\mathrm{tr}(\mathrm{Cov}(E)^k)\\ and \\h_0 = 1 - \frac{2 \theta_1
\theta_3}{3 \theta_2^2}\\.

## Usage

``` r
q_statistic_ucl(E, alpha = 0.01)
```

## Arguments

- E:

  Numeric matrix of training-set PCA residuals

- alpha:

  Significance level (default 0.01 = 99% confidence per Safanelli)

## Value

Numeric scalar UCL_Q
