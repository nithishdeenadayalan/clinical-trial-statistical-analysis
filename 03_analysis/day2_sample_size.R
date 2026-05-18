# -----------------------------------------------
# Clinical Trial Statistical Analysis
# Day 2 — Sample Size and Power Calculations
# Nithish Deenadayalan
# -----------------------------------------------

library(pwr)
library(survival)

data(lung, package = "survival")

# -----------------------------------------------
# STEP 1 — Sample size for survival analysis
# -----------------------------------------------

# We want to detect a difference in survival between
# male and female patients
# Assumptions:
# - 80% power (standard in clinical trials)
# - 5% significance level (two-sided)
# - Medium effect size (0.5 — standard assumption)

cat("--- Sample Size Calculation ---\n")
cat("Objective: Detect survival difference between male and female patients\n")
cat("Power: 80%\n")
cat("Significance level: 5% (two-sided)\n\n")

# Two sample comparison
result <- pwr.2p.test(
  h = 0.5,        # medium effect size
  sig.level = 0.05,
  power = 0.80,
  alternative = "two.sided"
)

print(result)
cat("\nMinimum sample size required per group:", ceiling(result$n), "\n")
cat("Total minimum sample size:", ceiling(result$n) * 2, "\n")
cat("Our dataset has:", nrow(lung), "patients — sufficient for analysis\n")

# -----------------------------------------------
# STEP 2 — Power of our actual study
# -----------------------------------------------

cat("\n--- Power of Our Actual Study ---\n")
cat("With", nrow(lung), "patients split by sex:\n")

male_n <- sum(lung$sex == 1, na.rm = TRUE)
female_n <- sum(lung$sex == 2, na.rm = TRUE)

cat("Males:", male_n, "\n")
cat("Females:", female_n, "\n")

# Calculate actual power with our sample size
actual_power <- pwr.2p2n.test(
  h = 0.5,
  n1 = male_n,
  n2 = female_n,
  sig.level = 0.05,
  alternative = "two.sided"
)

cat("Achieved power:", round(actual_power$power * 100, 1), "%\n")
cat("Conclusion: Study is", 
    ifelse(actual_power$power >= 0.80, "adequately powered", "underpowered"),
    "for detecting a medium effect size\n")

# -----------------------------------------------
# STEP 3 — Sensitivity analysis
# -----------------------------------------------

cat("\n--- Power at different effect sizes ---\n")
effect_sizes <- c(0.2, 0.3, 0.5, 0.8)
effect_labels <- c("Small", "Small-Medium", "Medium", "Large")

for(i in seq_along(effect_sizes)){
  p <- pwr.2p2n.test(
    h = effect_sizes[i],
    n1 = male_n,
    n2 = female_n,
    sig.level = 0.05,
    alternative = "two.sided"
  )
  cat(effect_labels[i], "effect (h=", effect_sizes[i], "):",
      round(p$power * 100, 1), "% power\n")
}