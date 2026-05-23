# Clinical Trial Statistical Analysis
### Survival Analysis of Lung Cancer Patients — NCCTG Clinical Trial
**Author:** Nithish Deenadayalan | MSc Business Analytics, UCD Smurfit Graduate Business School

[![Live Dashboard](https://img.shields.io/badge/Live%20Dashboard-Streamlit-red)](YOUR_STREAMLIT_LINK)
[![GitHub](https://img.shields.io/badge/GitHub-Repository-black)](https://github.com/nithishdeenadayalan/clinical-trial-statistical-analysis)

---

## Overview

This project demonstrates a complete end-to-end statistical analysis pipeline on a real clinical trial dataset, modelled on the analytical workflows used in pharmaceutical and clinical research organisations such as ICON plc, Parexel, and Quintiles.

The analysis covers every stage of a clinical trial statistical workflow — from data documentation and statistical analysis planning through to survival analysis, Cox regression, repeated measures modelling, and interactive dashboard delivery.

---

## Live Dashboard

The interactive dashboard is deployed publicly on Streamlit Cloud and allows users to filter patients by sex, ECOG performance score, and age range, with all analyses updating in real time.

**Live link:** YOUR_STREAMLIT_LINK

---

## Dataset

The NCCTG Lung Cancer dataset is a real clinical trial dataset from the North Central Cancer Treatment Group, available in R's survival package.

- 228 lung cancer patients
- 10 clinical variables including survival time, vital status, age, sex, and ECOG performance score
- 165 events (deaths), 63 censored observations
- Median survival: 310 days

---

## Statistical Methods

### 1. Descriptive Statistics
Baseline characteristics were summarised overall and stratified by sex, including means, standard deviations, frequency counts, and p-values using the gtsummary package in R.

### 2. Sample Size and Power Calculations
Minimum sample size and achieved power were calculated using the pwr package in R. The dataset achieves 95.8% power for detecting a medium effect size (Cohen's h = 0.5) at a two-sided 5% significance level.

### 3. Kaplan Meier Survival Analysis
Survival functions were estimated using the Kaplan Meier method for the overall cohort and stratified by sex and ECOG performance score. Group differences were assessed using the log-rank test.

### 4. Cox Proportional Hazards Regression
A multivariate Cox regression model was fitted with age, sex, and ECOG score as covariates. The proportional hazards assumption was verified using Schoenfeld residuals. Analyses were conducted in R and replicated in Python using the lifelines library.

### 5. Repeated Measures Analysis
A linear mixed effects model was fitted using lme4 in R to assess ECOG performance score trajectory over time, with patient included as a random effect to account for within-patient correlation.

---

## Key Findings

- Overall median survival was 310 days (95% CI: varies by subgroup)
- Female patients had a 42% lower hazard of death compared to male patients (HR 0.58, 95% CI 0.41 to 0.80, p < 0.001)
- ECOG performance score was the strongest independent predictor of survival — each one unit increase raised the hazard of death by 59% (HR 1.59, 95% CI 1.27 to 2.0, p < 0.001)
- Age was not a significant independent predictor after adjusting for sex and ECOG score (HR 1.01, p = 0.232)
- Sex difference in survival was statistically significant by log-rank test (p = 0.0016)
- ECOG group difference in survival was highly significant (p < 0.0001)

---

## Technologies Used

| Category | Tools |
|---|---|
| Statistical Analysis | R (survival, survminer, lme4, pwr, gtsummary) |
| Data Processing | Python (pandas, numpy, lifelines) |
| Visualisation | Plotly, ggplot2, survminer |
| Dashboard | Streamlit |
| SAS Replication | SAS OnDemand (PROC LIFETEST, PROC PHREG, PROC MIXED) |
| Version Control | Git, GitHub |

---

## How to Run Locally

**Clone the repository:**
```bash
git clone https://github.com/nithishdeenadayalan/clinical-trial-statistical-analysis.git
cd clinical-trial-statistical-analysis
```

**Install Python dependencies:**
```bash
pip install -r 05_dashboard/requirements.txt
```

**Run the dashboard:**
```bash
cd 05_dashboard
streamlit run app.py
```

**Run the R analysis scripts:**

Open RStudio, set your working directory to the project root, and run the scripts in the `03_analysis` folder in order from day1 through day5.

---

## About the Author

LinkedIn: linkedin.com/in/nithish-deenadayalan
Email: nitd020611@gmail.com
