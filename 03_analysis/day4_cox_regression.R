# -----------------------------------------------
# Clinical Trial Statistical Analysis
# Day 4 — Cox Proportional Hazards Regression
# Nithish Deenadayalan
# -----------------------------------------------

library(survival)
library(survminer)
library(dplyr)

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
  filter(!is.na(ph.ecog), !is.na(age), !is.na(sex))

surv_obj <- Surv(time = lung_clean$time, event = lung_clean$status == 2)

# -----------------------------------------------
# STEP 2 — Univariate Cox models
# -----------------------------------------------

cat("--- Univariate Cox Regression ---\n\n")

covariates <- c("age", "sex", "ph.ecog")
for(var in covariates){
  formula  <- as.formula(paste("surv_obj ~", var))
  fit      <- coxph(formula, data = lung_clean)
  s        <- summary(fit)
  hr       <- round(s$conf.int[1], 3)
  ci_lower <- round(s$conf.int[3], 3)
  ci_upper <- round(s$conf.int[4], 3)
  pval     <- round(s$coefficients[5], 4)
  cat(var, "— HR:", hr, "| 95% CI:", ci_lower, "to", ci_upper, "| p-value:", pval, "\n")
}

# -----------------------------------------------
# STEP 3 — Multivariate Cox model
# -----------------------------------------------

cat("\n--- Multivariate Cox Regression ---\n")

cox_multi <- coxph(surv_obj ~ age + sex + ph.ecog, data = lung_clean)
summary(cox_multi)

# -----------------------------------------------
# STEP 4 — Extract and display results table
# -----------------------------------------------

cat("\n--- Results Table ---\n")

cox_summary <- summary(cox_multi)

results <- data.frame(
  Covariate = c("Age (per year)", "Sex (Female vs Male)", "ECOG Score"),
  HR        = round(cox_summary$conf.int[,1], 3),
  CI_Lower  = round(cox_summary$conf.int[,3], 3),
  CI_Upper  = round(cox_summary$conf.int[,4], 3),
  P_Value   = round(cox_summary$coefficients[,5], 4)
)

print(results)

cat("\nInterpretations:\n")
cat("Age: For every 1 year increase in age, the hazard of death changes by",
    round((results$HR[1]-1)*100, 1), "%\n")
cat("Sex: Females have", round((1-results$HR[2])*100, 1),
    "% lower hazard of death compared to males\n")
cat("ECOG: Each 1 unit increase in ECOG score increases hazard of death by",
    round((results$HR[3]-1)*100, 1), "%\n")

# -----------------------------------------------
# STEP 5 — Test proportional hazards assumption
# -----------------------------------------------

cat("\n--- Proportional Hazards Assumption Test ---\n")
ph_test <- cox.zph(cox_multi)
print(ph_test)
cat("\nInterpretation: p > 0.05 means assumption is not violated\n")
cat("Global p-value:", round(ph_test$table["GLOBAL","p"], 4), "\n")
cat("Assumption",
    ifelse(ph_test$table["GLOBAL","p"] > 0.05, "NOT violated — model is valid", 
           "violated — consider time-varying coefficients"), "\n")

# -----------------------------------------------
# STEP 6 — Forest plot of hazard ratios
# -----------------------------------------------

p_forest <- ggforest(
  cox_multi,
  data      = lung_clean,
  main      = "Hazard Ratios — Multivariate Cox Regression",
  cpositions = c(0.02, 0.22, 0.4),
  fontsize  = 0.9
)
print(p_forest)

# -----------------------------------------------
# STEP 7 — Schoenfeld residuals plot
# -----------------------------------------------

png("03_analysis/plot_schoenfeld_residuals.png", width=900, height=600, res=120)
ggcoxzph(ph_test)
dev.off()

# -----------------------------------------------
# STEP 8 — Save results
# -----------------------------------------------

png("03_analysis/plot_forest_cox.png", width=900, height=500, res=120)
ggforest(cox_multi, data = lung_clean,
         main = "Hazard Ratios — Multivariate Cox Regression")
dev.off()

write.csv(results, "03_analysis/cox_regression_results.csv", row.names = FALSE)

cat("\nForest plot and results saved to 03_analysis folder\n")
cat("\nDay 4 complete. Ready for Day 5 — Repeated Measures Analysis.\n")