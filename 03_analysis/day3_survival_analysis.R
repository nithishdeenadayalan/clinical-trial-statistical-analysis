# -----------------------------------------------
# Clinical Trial Statistical Analysis
# Day 3 — Kaplan Meier Survival Analysis
# Nithish Deenadayalan
# -----------------------------------------------

library(survival)
library(survminer)
library(dplyr)
library(ggplot2)

data(lung, package = "survival")

# -----------------------------------------------
# STEP 1 — Prepare data
# -----------------------------------------------

lung_clean <- lung %>%
  mutate(
    sex_label  = factor(sex, levels = c(1,2), labels = c("Male","Female")),
    ecog_label = factor(ph.ecog, levels = c(0,1,2,3),
                        labels = c("Fully Active","Restricted","Ambulatory","Limited"))
  ) %>%
  filter(!is.na(ph.ecog))

# -----------------------------------------------
# STEP 2 — Overall survival object
# -----------------------------------------------

surv_obj <- Surv(time = lung_clean$time, event = lung_clean$status == 2)

cat("--- Overall Survival Summary ---\n")
overall_fit <- survfit(surv_obj ~ 1, data = lung_clean)
print(summary(overall_fit)$table)

cat("\nMedian survival:", summary(overall_fit)$table["median"], "days\n")
cat("95% CI:", summary(overall_fit)$table["0.95LCL"], 
    "to", summary(overall_fit)$table["0.95UCL"], "days\n")

# -----------------------------------------------
# STEP 3 — Survival by sex
# -----------------------------------------------

cat("\n--- Survival by Sex ---\n")
fit_sex <- survfit(surv_obj ~ sex_label, data = lung_clean)
print(summary(fit_sex)$table)

# Log rank test
logrank_sex <- survdiff(surv_obj ~ sex_label, data = lung_clean)
cat("\nLog-rank test p-value:", 
    round(1 - pchisq(logrank_sex$chisq, df=1), 4), "\n")
cat("Interpretation:", 
    ifelse(1 - pchisq(logrank_sex$chisq, df=1) < 0.05,
           "Significant survival difference between sexes",
           "No significant survival difference between sexes"), "\n")

# -----------------------------------------------
# STEP 4 — Survival by ECOG score
# -----------------------------------------------

cat("\n--- Survival by ECOG Score ---\n")
fit_ecog <- survfit(surv_obj ~ ecog_label, data = lung_clean)
print(summary(fit_ecog)$table)

logrank_ecog <- survdiff(surv_obj ~ ecog_label, data = lung_clean)
cat("\nLog-rank test p-value:", 
    round(1 - pchisq(logrank_ecog$chisq, df=3), 6), "\n")

# -----------------------------------------------
# STEP 5 — Plot overall survival curve
# -----------------------------------------------

p1 <- ggsurvplot(
  overall_fit,
  data          = lung_clean,
  conf.int      = TRUE,
  risk.table    = TRUE,
  surv.median.line = "hv",
  title         = "Overall Survival — NCCTG Lung Cancer Trial",
  xlab          = "Time (Days)",
  ylab          = "Survival Probability",
  ggtheme       = theme_minimal(),
  palette       = "#2E75B6",
  risk.table.height = 0.25,
  legend        = "none"
)
print(p1)

# -----------------------------------------------
# STEP 6 — Plot survival by sex
# -----------------------------------------------

p2 <- ggsurvplot(
  fit_sex,
  data             = lung_clean,
  conf.int         = TRUE,
  risk.table       = TRUE,
  pval             = TRUE,
  pval.method      = TRUE,
  surv.median.line = "hv",
  title            = "Survival by Sex — NCCTG Lung Cancer Trial",
  xlab             = "Time (Days)",
  ylab             = "Survival Probability",
  ggtheme          = theme_minimal(),
  palette          = c("#2E75B6","#E74C3C"),
  legend.title     = "Sex",
  legend.labs      = c("Male","Female"),
  risk.table.height = 0.25
)
print(p2)

# -----------------------------------------------
# STEP 7 — Plot survival by ECOG score
# -----------------------------------------------

p3 <- ggsurvplot(
  fit_ecog,
  data             = lung_clean,
  conf.int         = FALSE,
  risk.table       = TRUE,
  pval             = TRUE,
  pval.method      = TRUE,
  surv.median.line = "hv",
  title            = "Survival by ECOG Performance Score",
  xlab             = "Time (Days)",
  ylab             = "Survival Probability",
  ggtheme          = theme_minimal(),
  palette          = c("#27AE60","#2E75B6","#E74C3C","#8E44AD"),
  legend.title     = "ECOG Score",
  legend.labs      = c("Fully Active","Restricted","Ambulatory","Limited"),
  risk.table.height = 0.30
)
print(p3)

# -----------------------------------------------
# STEP 8 — Save plots
# -----------------------------------------------

# Save each plot as PNG
png("03_analysis/plot_overall_survival.png", width=900, height=700, res=120)
print(p1)
dev.off()

png("03_analysis/plot_survival_by_sex.png", width=900, height=700, res=120)
print(p2)
dev.off()

png("03_analysis/plot_survival_by_ecog.png", width=900, height=700, res=120)
print(p3)
dev.off()

cat("\nAll 3 survival plots saved to 03_analysis folder\n")
cat("\nKey findings:\n")
cat("- Overall median survival:", summary(overall_fit)$table["median"], "days\n")
cat("- Log-rank sex p-value:", round(1 - pchisq(logrank_sex$chisq, df=1), 4), "\n")
cat("- Log-rank ECOG p-value:", round(1 - pchisq(logrank_ecog$chisq, df=3), 6), "\n")
cat("\nDay 3 complete. Ready for Day 4 — Cox Regression.\n")