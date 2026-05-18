# -----------------------------------------------
# Clinical Trial Statistical Analysis
# Day 2 — Descriptive Statistics
# Nithish Deenadayalan
# -----------------------------------------------

library(survival)
library(dplyr)
library(gtsummary)

data(lung, package = "survival")

# -----------------------------------------------
# STEP 1 — Clean the data
# -----------------------------------------------

lung_clean <- lung %>%
  mutate(
    sex_label  = factor(sex,      levels = c(1,2), labels = c("Male","Female")),
    ecog_label = factor(ph.ecog,  levels = c(0,1,2,3),
                        labels = c("Fully Active","Restricted","Ambulatory","Limited")),
    age_group  = cut(age, breaks = c(0,55,65,75,100),
                     labels = c("Under 55","55-64","65-74","75+"),
                     right  = FALSE),
    status_label = factor(status, levels = c(1,2),
                          labels = c("Censored","Dead"))
  )

cat("Data cleaned successfully\n")
cat("Rows:", nrow(lung_clean), "\n")

# -----------------------------------------------
# STEP 2 — Overall summary
# -----------------------------------------------

cat("\n--- Overall Summary ---\n")
cat("Total patients:      ", nrow(lung_clean), "\n")
cat("Deaths (events):     ", sum(lung_clean$status == 2), "\n")
cat("Censored:            ", sum(lung_clean$status == 1), "\n")
cat("Event rate:          ", round(mean(lung_clean$status == 2)*100,1), "%\n")
cat("Median survival:     ", median(lung_clean$time), "days\n")
cat("Mean age:            ", round(mean(lung_clean$age),1), "years\n")
cat("Age range:           ", min(lung_clean$age), "to", max(lung_clean$age), "years\n")

# -----------------------------------------------
# STEP 3 — Summary by sex
# -----------------------------------------------

cat("\n--- Summary by Sex ---\n")
sex_summary <- lung_clean %>%
  group_by(sex_label) %>%
  summarise(
    n              = n(),
    deaths         = sum(status == 2),
    event_rate_pct = round(mean(status == 2)*100, 1),
    median_age     = round(median(age), 1),
    median_survival= median(time),
    mean_survival  = round(mean(time), 1)
  )
print(sex_summary)

# -----------------------------------------------
# STEP 4 — Summary by ECOG score
# -----------------------------------------------

cat("\n--- Summary by ECOG Performance Score ---\n")
ecog_summary <- lung_clean %>%
  filter(!is.na(ecog_label)) %>%
  group_by(ecog_label) %>%
  summarise(
    n               = n(),
    deaths          = sum(status == 2),
    event_rate_pct  = round(mean(status == 2)*100, 1),
    median_survival = median(time),
    mean_survival   = round(mean(time), 1)
  )
print(ecog_summary)

# -----------------------------------------------
# STEP 5 — Age group breakdown
# -----------------------------------------------

cat("\n--- Age Group Distribution ---\n")
age_summary <- lung_clean %>%
  group_by(age_group) %>%
  summarise(
    n               = n(),
    pct             = round(n()/nrow(lung_clean)*100, 1),
    deaths          = sum(status == 2),
    median_survival = median(time)
  )
print(age_summary)

# -----------------------------------------------
# STEP 6 — Baseline characteristics table
# -----------------------------------------------

cat("\n--- Building Baseline Characteristics Table ---\n")

tbl <- lung_clean %>%
  select(age, sex_label, ecog_label, age_group,
         ph.karno, wt.loss, time, status_label) %>%
  tbl_summary(
    by = sex_label,
    label = list(
      age          ~ "Age (years)",
      ecog_label   ~ "ECOG Performance Score",
      age_group    ~ "Age Group",
      ph.karno     ~ "Karnofsky Score (Physician)",
      wt.loss      ~ "Weight Loss (lbs)",
      time         ~ "Survival Time (days)",
      status_label ~ "Vital Status"
    ),
    statistic = list(
      all_continuous()  ~ "{mean} ({sd})",
      all_categorical() ~ "{n} ({p}%)"
    ),
    missing = "ifany"
  ) %>%
  add_overall() %>%
  add_p() %>%
  bold_labels()

print(tbl)

# Save summary tables to CSV
write.csv(sex_summary,  "01_data/summary_by_sex.csv",  row.names = FALSE)
write.csv(ecog_summary, "01_data/summary_by_ecog.csv", row.names = FALSE)
write.csv(age_summary,  "01_data/summary_by_age.csv",  row.names = FALSE)

cat("\nAll summary tables saved to 01_data folder\n")
cat("\nDay 2 complete. Ready for Day 3 — Survival Analysis.\n")