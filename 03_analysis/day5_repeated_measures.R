# -----------------------------------------------
# Clinical Trial Statistical Analysis
# Day 5 — Repeated Measures Analysis
# Nithish Deenadayalan
# -----------------------------------------------

library(survival)
library(lme4)
library(lmerTest)
library(dplyr)
library(ggplot2)

data(lung, package = "survival")

# -----------------------------------------------
# STEP 1 — Prepare longitudinal data
# -----------------------------------------------

# We will simulate repeated ECOG measurements over time
# based on the available data to demonstrate repeated
# measures methodology as would be done in a real trial

set.seed(42)

lung_clean <- lung %>%
  filter(!is.na(ph.ecog), !is.na(age), !is.na(sex)) %>%
  mutate(
    patient_id = row_number(),
    sex_label  = factor(sex, levels = c(1,2), labels = c("Male","Female")),
    baseline_ecog = ph.ecog
  )

# Create repeated measurements at 3 time points
# Time 0 = baseline, Time 1 = 3 months, Time 2 = 6 months
long_data <- lung_clean %>%
  slice(rep(1:n(), each = 3)) %>%
  mutate(
    timepoint = rep(c(0, 90, 180), times = nrow(lung_clean)),
    ecog_score = baseline_ecog +
      (timepoint / 90) * 0.3 +
      rnorm(n(), mean = 0, sd = 0.3),
    ecog_score = pmax(0, pmin(3, round(ecog_score, 1)))
  )

cat("Longitudinal dataset created\n")
cat("Patients:", n_distinct(long_data$patient_id), "\n")
cat("Total observations:", nrow(long_data), "\n")
cat("Time points: 0, 90, 180 days\n")

# -----------------------------------------------
# STEP 2 — Descriptive summary by timepoint
# -----------------------------------------------

cat("\n--- ECOG Score by Timepoint ---\n")
time_summary <- long_data %>%
  group_by(timepoint) %>%
  summarise(
    n           = n(),
    mean_ecog   = round(mean(ecog_score), 3),
    sd_ecog     = round(sd(ecog_score), 3),
    median_ecog = round(median(ecog_score), 3)
  )
print(time_summary)

cat("\n--- ECOG Score by Timepoint and Sex ---\n")
sex_time_summary <- long_data %>%
  group_by(timepoint, sex_label) %>%
  summarise(
    n         = n(),
    mean_ecog = round(mean(ecog_score), 3),
    sd_ecog   = round(sd(ecog_score), 3),
    .groups   = "drop"
  )
print(sex_time_summary)

# -----------------------------------------------
# STEP 3 — Linear mixed effects model
# -----------------------------------------------

cat("\n--- Linear Mixed Effects Model ---\n")
cat("Model: ECOG score ~ time + sex + age + (1|patient)\n\n")

lme_model <- lmer(
  ecog_score ~ timepoint + sex_label + age + (1 | patient_id),
  data = long_data,
  REML = TRUE
)

print(summary(lme_model))

# -----------------------------------------------
# STEP 4 — Extract fixed effects
# -----------------------------------------------

cat("\n--- Fixed Effects Summary ---\n")
fixed_effects <- as.data.frame(coef(summary(lme_model)))
fixed_effects$term <- rownames(fixed_effects)
fixed_effects <- fixed_effects %>%
  select(term, Estimate, `Std. Error`, `t value`, `Pr(>|t|)`) %>%
  mutate(across(where(is.numeric), ~ round(., 4)))

print(fixed_effects)

cat("\nInterpretation:\n")
time_coef <- fixed_effects$Estimate[fixed_effects$term == "timepoint"]
sex_coef  <- fixed_effects$Estimate[fixed_effects$term == "sex_labelFemale"]
cat("Time effect: ECOG score changes by", round(time_coef * 90, 4),
    "units per 90 days\n")
cat("Sex effect: Females have", round(sex_coef, 4),
    "lower ECOG score than males on average\n")

# -----------------------------------------------
# STEP 5 — Random effects
# -----------------------------------------------

cat("\n--- Random Effects (Between-patient variability) ---\n")
re_var <- as.data.frame(VarCorr(lme_model))
print(re_var)
cat("ICC (Intraclass Correlation):",
    round(re_var$vcov[1] / sum(re_var$vcov), 3), "\n")
cat("Interpretation: proportion of total variance due to between-patient differences\n")

# -----------------------------------------------
# STEP 6 — Plot ECOG trajectories over time
# -----------------------------------------------

p1 <- ggplot(time_summary, aes(x = timepoint, y = mean_ecog)) +
  geom_line(color = "#2E75B6", linewidth = 1.2) +
  geom_point(color = "#2E75B6", size = 4) +
  geom_errorbar(aes(ymin = mean_ecog - sd_ecog,
                    ymax = mean_ecog + sd_ecog),
                width = 10, color = "#2E75B6", alpha = 0.6) +
  labs(
    title   = "Mean ECOG Performance Score Over Time",
    x       = "Time (Days)",
    y       = "Mean ECOG Score",
    caption = "Error bars represent ±1 SD"
  ) +
  scale_x_continuous(breaks = c(0, 90, 180)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p1)

# -----------------------------------------------
# STEP 7 — Plot by sex
# -----------------------------------------------

p2 <- ggplot(sex_time_summary,
             aes(x = timepoint, y = mean_ecog,
                 color = sex_label, group = sex_label)) +
  geom_line(linewidth = 1.2) +
  geom_point(size = 4) +
  geom_errorbar(aes(ymin = mean_ecog - sd_ecog,
                    ymax = mean_ecog + sd_ecog),
                width = 10, alpha = 0.6) +
  scale_color_manual(values = c("#2E75B6","#E74C3C")) +
  labs(
    title   = "ECOG Performance Score Over Time by Sex",
    x       = "Time (Days)",
    y       = "Mean ECOG Score",
    color   = "Sex",
    caption = "Error bars represent ±1 SD"
  ) +
  scale_x_continuous(breaks = c(0, 90, 180)) +
  theme_minimal() +
  theme(plot.title = element_text(size = 14, face = "bold"))

print(p2)

# -----------------------------------------------
# STEP 8 — Save outputs
# -----------------------------------------------

png("03_analysis/plot_ecog_trajectory.png", width=900, height=600, res=120)
print(p1)
dev.off()

png("03_analysis/plot_ecog_by_sex.png", width=900, height=600, res=120)
print(p2)
dev.off()

write.csv(time_summary,     "03_analysis/repeated_measures_summary.csv",    row.names=FALSE)
write.csv(sex_time_summary, "03_analysis/repeated_measures_by_sex.csv",     row.names=FALSE)
write.csv(fixed_effects,    "03_analysis/repeated_measures_fixed_effects.csv", row.names=FALSE)

cat("\nAll repeated measures outputs saved\n")
cat("\nDay 5 complete. Ready for Day 6 — Streamlit Dashboard.\n")