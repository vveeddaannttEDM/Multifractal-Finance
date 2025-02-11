# Multifractal-Finance
# Multifractal Analysis of Financial Volatility

This project implements the multifractal analysis of financial time series inspired by the paper *Direct Evidence for Inversion Formula in Multifractal Financial Volatility Measure* by Zhi-Qiang Jiang and Wei-Xing Zhou. The goal is to verify the inversion formula in financial markets using high-frequency S&P 500 data.

Rolling Window Analysis:
This section breaks the volatility series into overlapping time windows and estimates the scaling exponent 
q=2) in each window. This shows how the multifractal spectrum might evolve over time.

Surrogate Analysis:
By randomly shuffling (surrogating) the volatility data, we “destroy” any temporal correlations. Recomputing the scaling exponents on the surrogate data tests whether the observed multifractality is due to genuine correlations in the data.

Bootstrap Resampling for Confidence Intervals:
A bootstrap procedure is implemented to obtain 95% confidence intervals for the scaling exponents 
. This gives a measure of the statistical uncertainty in the estimation.
## Setup

1. Clone the repository:
```bash
git clone https://github.com/your-username/multifractal_finance_analysis.git
cd multifractal_finance_analysis


