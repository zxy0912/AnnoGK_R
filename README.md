# GhostKnockoff Pipeline Example

This example demonstrates how to use the `AnnoGKR` package to perform variable selection on simulated GWAS/TWAS data, comparing the performance of standard GhostKnockoff against annotation-informed versions.

## 1. Data Simulation
We simulate a linear trait with $N=5000$ samples and $p=600$ features, assuming an AR(1) correlation structure for the features.

```R
library(AnnoGKR)

# Simulation Parameters
seed <- 1234
set.seed(seed)
N.effect <- 5000; p <- 600; n0 <- 30; heri <- 0.1; rho <- 0.5
M <- 1; threshold <- 0.1
amplitude <- sqrt(heri / (1 - heri))

# True signal generation
sigprob <- rep(0, p); sigprob[1:60] <- 1/(1:60)^2/sum(1/(1:60)^2)
rand <- sample(1:p, n0, prob = sigprob)
rand_sign <- sample(c(-1, 1), size = p, replace = TRUE)

# Covariates and Response
Covariance <- toeplitz(rho^(0:(p-1)))
X <- scale(matrix(rnorm(N.effect*p), N.effect, p) %*% chol(Covariance))
beta <- rep(0, p);
beta[rand] = amplitude
beta = beta*rand_sign
y <- scale(X %*% beta + sqrt(n0) * rnorm(N.effect))

# Summary Statistics and LD Matrix
Z <- sapply(1:p, function(l) coef(summary(lm(y ~ X[,l])))[2,3])
LD <- cor(X)
R  <- scale(as.matrix(1:p)) # Annotation matrix
```
## 2. Comparison of Methods
We evaluate three strategies to see how incorporating structural annotations improves statistical power while controlling FDR.

Method A: Standard GhostKnockoff (No Annotations)

```R
set.seed(seed)
fit.prelim <- GhostKnockoff.prelim(cor.G = LD, M = M, method = "sdp")
GK1_lasso  <- GhostKnockoff.fit(Z, N.effect, fit.prelim, method = 'lasso')
GK.filter  <- GhostKnockoff.filter(GK1_lasso$T_0[[1]], GK1_lasso$T_k[[1]])

rej <- which(GK.filter$q <= threshold)
cat("Standard Power:", power_cal(rej, rand), "FDR:", fdr_cal(rej, rand), "\n")
```

Method B: AnnoGK-Simple (Annotation-Informed)

```R
set.seed(seed)
GK_simple_res <- GK_simple(Z=Z, R=R, M=M, LD=LD, n=N.effect, ts='lasso')
beta <- GK_simple_res$beta
GK.filter <- GhostKnockoff.filter(abs(beta[1:p]), abs(matrix(beta[-(1:p)], p, M)))

rej <- which(GK.filter$q <= threshold)
cat("Simple-Anno Power:", power_cal(rej, rand), "FDR:", fdr_cal(rej, rand), "\n")
```

Method C: Full AnnoGK-M (Multi-copy Annotation)

```R
set.seed(seed)
GK_anno_res <- GK_anno_M(Z, R, M, LD, N.effect)
GK.filter   <- GhostKnockoff.filter(GK_anno_res$T_0, GK_anno_res$T_k)

rej <- which(GK.filter$q <= threshold)
cat("Full-Anno Power:", power_cal(rej, rand), "FDR:", fdr_cal(rej, rand), "\n")
```

## 3. Performance Comparison

We evaluate three strategies to see how incorporating structural annotations improves statistical power while controlling FDR.

### Performance Summary

| Method | Power | FDR |
| :--- | :--- | :--- |
| **Standard GhostKnockoff** | 0.000 | 0.000 |
| **AnnoGK-Simple** | **0.867** | **0.000** |
| **Full AnnoGK-M** | 0.767 | 0.000 |
