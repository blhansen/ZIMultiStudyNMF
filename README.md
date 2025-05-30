# ZIMultiStudyNMF
Methods for Sampling from the Zero Inflated Bayesian Multi-Study Non-negative Matrix Factorization model

# Installation
Installation can be done via devtools:
```{r installation, echo = TRUE, results = TRUE, tidy = TRUE}
devtools::install_github('blhansen/ZIMultiStudyNMF')
library(ZIMSNMF)
```

# Example with simulated data:
```{r example, echo = TRUE, results = TRUE, tidy = TRUE}
S <- 3
P <- 20
K <- 5
N_s <- c(50, 60, 70)

W <- matrix(rgamma(P*K, shape=1, scale=2), nrow=P, ncol=K)
H_s <- lapply(1:S, function(s){matrix(rgamma(K*N_s[s], shape=10, rate=10/20), nrow=K, ncol=N_s[s])})
M_s <- lapply(1:S, function(s) matrix(rpois(N_s[s]*P, W %*% H_s[[s]], nrow=K, ncol=N_s[s]))

gibbs_zinmf(M_s=M_s, K=5)
```
